from textwrap import dedent
from LookupTable import LUT
import paths, os, sys, subprocess, datetime
import json
import pdb

class PatternOptimization:
    def __init__(self, config):
        self.dim = config['dim']
        self.pattern = config['pattern']
        self.material = paths.material(config['material'])
        self.patoptArgs = config.get('args', ['-I', '--solver', 'slsqp',
                                              '--TensorFitConstraint',
                                              '--tensor_fit_tolerance=1e-5', '-n50'])

        self.jobTemplate = None
        if 'jobTemplate' in config:
            jtpath = config['jobTemplate']
            f = open(jtpath, 'r')
            self.jobTemplate = json.load(open(config['jobTemplate']))

        # Read variable bounds from config, with dim-dependent defaults
        if (self.dim == 2):
            self.radiusBounds      = config.get(     'radiusBounds', [0.02, 0.17])
            self.translationBounds = config.get('translationBounds', [-0.14, 0.14])
            self.blendingBounds    = config.get('blendingBounds', [0.0078125, 0.2])
        elif (self.dim == 3):
            self.radiusBounds      = config.get(     'radiusBounds', [0.02, 0.2])
            self.translationBounds = config.get('translationBounds', [-0.2, 0.2])
            self.blendingBounds    = config.get('blendingBounds', [0.0078125, 0.2])
        else: raise Exception("Invalid dimension")

        self.jobs = []
        self.outputs = []
        self.constraints = []

    def writeJobs(self, directory):
        os.mkdir(directory)
        for i, j in enumerate(self.jobs):
            f = open(directory + '/%i.job' % i, 'w')
            f.write(j)
            f.close()

    # expects a list of constraints in string form, e.g.:
    # ["p3 = p5 + 1.0 / 6.0", "p4 = p7 + 10 / 6.0"]
    def setConstraints(self, constraints):
        self.constraints = constraints

    # Create and enqueue a job for a particular target tensor, initialParams
    # Use template json object self.jobTemplate if it exists
    def enqueueJob(self, targetE, targetNu, initialParams):
        # TODO: make bounds configurable

        formatSequence = lambda x: ", ".join(map(str, x))

        constraintString = ""
        if len(self.constraints) > 0:
            constraintString = '\n"paramConstraints": [%s],' % ", ".join(['"%s"' % c for c in self.constraints])
        jobContents = ""
        if self.jobTemplate != None:
            self.jobTemplate['target'] = {'type': 'isotropic', 'young': targetE, 'poisson': targetNu};
            self.jobTemplate['initial_params'] = list(initialParams);
            jobContents = json.dumps(self.jobTemplate, indent=4)
        else:
            jobContents = """\
            {
                "dim": %i,
                "target": {
                    "type": "isotropic",
                    "young": %f,
                    "poisson": %f
                },%s
                "initial_params": [%s],
                "radiusBounds": [%s],
                "translationBounds": [%s],
                "blendingBounds": [%s]
            }
            """
            jobContents = dedent(jobContents) % (self.dim, targetE, targetNu, constraintString,
                    formatSequence(initialParams), formatSequence(self.radiusBounds),
                    formatSequence(self.translationBounds),
                    formatSequence(self.blendingBounds))
        self.jobs.append(jobContents)
   

    def run(self, directory):
        self.writeJobs(directory)
        for i in range(len(self.jobs)):
            itPrefix = directory + '/job-%i_it' % i
            cmd = [paths.optimizer(self.dim), directory + '/%i.job' % i,
                '-p', paths.pattern(self.pattern, self.dim), '-m', self.material, '-o', itPrefix] + self.patoptArgs
            outPath = directory + '/stdout_%i.txt' % i
            with open(outPath, 'w') as outLog:
                try:
                    ret = subprocess.call(cmd, stdout=outLog)
                    if (ret != 0):
                        print cmd
                        print "FAILED optimization (nonzero exit status)"
                except:
                    print cmd
                    sys.stderr.write("WARNING: optimization '%s/%i.job' died\n" % (directory, i))
        self.jobs = [] # remove finished jobs from queue

    # Submit an array job for jobs numbered firstJobIndex..lastJobIndex
    # Return the array job's id.
    def submitArrayJob(self, directory, firstJobIndex, lastJobIndex, nprocs, walltime, mem):
        if (lastJobIndex < firstJobIndex): raise Exception("Invalid job index range")
        cmd = [paths.optimizer(self.dim), directory + '/${PBS_ARRAYID}.job',
            '-p', paths.pattern(self.pattern, self.dim), '-m', self.material] + self.patoptArgs
        # We need to kill our job a couple minutes before the walltime expires so
        # that we can recover stdout/stderr. Subtract 2 minutes from the requested walltime.
        # parse H:M:S from walltime (python's datetime/time classes are annoyingly incomplete)
        hms = walltime.split(':')
        e = Exception("Error parsing walltime; must be in the format h:m:s")
        if (len(hms) != 3): raise e
        h,m,s = map(int, hms)
        if (h > 24 | h < 0 | m > 60 | m < 0 | s > 60 | s < 0): raise e
        fullTime = datetime.timedelta(hours=h, minutes=m, seconds=s)
        killTime = fullTime - datetime.timedelta(minutes=2)

        # TODO: figure out how to redirect stdout/sterr out of existence

        pbsScript = """\
        ###-----PBS Directives Start-----###

        #PBS -V
        #PBS -S /bin/bash
        #PBS -N {name}
        #PBS -l nodes=1:ppn={nprocs}
        #PBS -l walltime={walltime}
        #PBS -l mem={mem}
        #PBS -M fjp234@nyu.edu
        #PBS -m a
        #PBS -j oe
        #PBS -o /state/partition1/${{PBS_JOBNAME}}.o${{PBS_JOBID}}.${{PBS_ARRAYID}}
        #PBS -t {firstIndex}-{lastIndex}

        ###-----PBS Directives End-----###
        cd ${{PBS_O_WORKDIR}}

        STDOUT_FILE=${{PBS_MEMDISK}}/stdout_${{PBS_JOBID}}.${{PBS_ARRAYID}}.txt
        timeout -s KILL {killseconds} {command} > $STDOUT_FILE 2>&1

        mv $STDOUT_FILE {directory}/stdout_${{PBS_ARRAYID}}.txt
        """
        pbsScript = dedent(pbsScript).format(
                name="ac_%s_%s" % (str(self.pattern), directory),
                firstIndex=firstJobIndex, lastIndex=lastJobIndex,
                command = " ".join(cmd), directory=directory,
                nprocs=nprocs, walltime=walltime, mem=mem,
                killseconds=killTime.seconds)
        p = subprocess.Popen(["qsub", "-"], stdin=subprocess.PIPE,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
        qsub_stdout, qsub_stderr = p.communicate(pbsScript)
        if qsub_stderr != "":
            print "WARNING: nonempty stderr response from qsub: ", qsub_stderr
        return qsub_stdout
