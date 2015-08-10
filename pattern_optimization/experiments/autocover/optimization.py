from textwrap import dedent
from LookupTable import LUT
import paths, os, sys, subprocess

class PatternOptimization:
    def __init__(self, config):
        self.dim = config['dim']
        self.pattern = config['pattern']
        self.material = paths.material(config['material'])
        self.patoptArgs = config.get('args', ['-I', '--solver', 'levenberg_marquardt'])

        # Read variable bounds from config, with dim-dependent defaults
        if (self.dim == 2):
            self.radiusBounds      = config.get(     'radiusBounds', [0.02, 0.17])
            self.translationBounds = config.get('translationBounds', [-0.14, 0.14])
        elif (self.dim == 3):
            self.radiusBounds      = config.get(     'radiusBounds', [0.3, 0.8])
            self.translationBounds = config.get('translationBounds', [-0.3, 0.3])
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

    def enqueueJob(self, targetE, targetNu, initialParams):
        # TODO: make bounds configurable

        formatSequence = lambda x: ", ".join(map(str, x))

        constraintString = ""
        if len(self.constraints) > 0:
            constraintString = '\n"paramConstraints": [%s],' % ", ".join(['"%s"' % c for c in self.constraints])
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
            "translationBounds": [%s]
        }
        """
        jobContents = dedent(jobContents) % (self.dim, targetE, targetNu, constraintString,
                formatSequence(initialParams), formatSequence(self.radiusBounds),
                formatSequence(self.translationBounds))
        self.jobs.append(jobContents)

    def run(self, directory):
        self.writeJobs(directory)
        for i in range(len(self.jobs)):
            cmd = [paths.optimizer(self.dim), directory + '/%i.job' % i,
                '-p', paths.pattern(self.pattern, self.dim), '-m', self.material] + self.patoptArgs
            outPath = directory + '/stdout_%i.txt' % i
            with open(outPath, 'w') as outLog:
                try: subprocess.call(cmd, stdout=outLog)
                except: sys.stderr.write("WARNING: optimization '%s/%i.job' died\n" % (directory, i))
        self.jobs = [] # remove finished jobs from queue

    # Submit an array job for jobs numbered firstJobIndex..lastJobIndex
    # Return the array job's id.
    def submitArrayJob(self, directory, firstJobIndex, lastJobIndex):
        if (lastJobIndex <= firstJobIndex): raise Exception("Invalid job index range")
        cmd = [paths.optimizer(self.dim), directory + '/${PBS_ARRAYID}.job',
            '-p', paths.pattern(self.pattern, self.dim), '-m', self.material] + self.patoptArgs
        # TODO: figure out how to redirect stdout/sterr out of existence
        pbsScript = """\
        ###-----PBS Directives Start-----###

        #PBS -V
        #PBS -S /bin/bash
        #PBS -N {name}
        #PBS -l nodes=1:ppn=2
        #PBS -l walltime=1:30:00
        #PBS -l mem=4GB
        #PBS -M fjp234@nyu.edu
        #PBS -m a
        #PBS -t {firstIndex}-{lastIndex}

        ###-----PBS Directives End-----###
        cd ${{PBS_O_WORKDIR}}

        STDOUT_FILE=${{PBS_MEMDISK}}/stdout_${{PBS_JOBID}}.${{PBS_ARRAYID}}.txt
        timeout -s KILL 85m {command} > $STDOUT_FILE 2>&1

        mv $STDOUT_FILE {directory}/stdout_${{PBS_ARRAYID}}.txt
        """
        pbsScript = dedent(pbsScript).format(
                name="ac_%i_%s" % (self.pattern, directory),
                firstIndex=firstJobIndex, lastIndex=lastJobIndex,
                command = " ".join(cmd), directory=directory)
        p = subprocess.Popen(["qsub", "-"], stdin=subprocess.PIPE,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
        qsub_stdout, qsub_stderr = p.communicate(pbsScript)
        if qsub_stderr != "":
            print "WARNING: nonempty stderr response from qsub: ", qsub_stderr
        return qsub_stdout
