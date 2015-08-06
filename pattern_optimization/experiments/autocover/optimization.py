from textwrap import dedent
from LookupTable import LUT
import paths, os, sys, subprocess

class PatternOptimization:
    def __init__(self, dim, pat, mat):
        self.dim = dim
        self.pattern = pat
        self.material = mat
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

        if (self.dim == 2):
            radiusBounds = [0.02, 0.17]
            translationBounds = [-0.14, 0.14]
        elif (self.dim == 3):
            radiusBounds = [0.3, 0.8]
            transRange   = [-0.3, 0.3]
        else: raise Exception("Invalid dimension")

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
                formatSequence(initialParams), formatSequence(radiusBounds),
                formatSequence(translationBounds))
        self.jobs.append(jobContents)

    def run(self, directory, args=['-I', '--solver', 'levenberg_marquardt']):
        self.writeJobs(directory)
        for i in range(len(self.jobs)):
            cmd = [paths.optimizer(self.dim), directory + '/%i.job' % i,
                '-p', paths.pattern(self.pattern, self.dim), '-m', self.material] + args
            outPath = directory + '/stdout_%i.txt' % i
            with open(outPath, 'w') as outLog:
                try: subprocess.call(cmd, stdout=outLog)
                except: sys.stderr.write("WARNING: optimization '%s/%i.job' died\n" % (directory, i))
        self.jobs = [] # remove finished jobs from queue

    # Submit all jobs as an array job to the HPC batch system.
    # Return the array job's id.
    def submitArrayJob(self, directory, args=['--solver', 'levenberg_marquardt']):
        numJobs = len(self.jobs)
        self.writeJobs(directory)
        cmd = [paths.optimizer(self.dim), '${PBS_ARRAYID}',
            '-p', paths.pattern(self.pattern, self.dim), '-m', self.material] + args
        outPath = directory + '/stdout_%i.txt' % i
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
        #PBS -e localhost:${{PBS_O_WORKDIR}}/{name}.e${{PBS_JOBID}}
        #PBS -o localhost:${{PBS_O_WORKDIR}}/{name}.o${{PBS_JOBID}}
        #PBS -t 0-{maxIdx}

        ###-----PBS Directives End-----###
        cd ${{PBS_O_WORKDIR}}
        {command}
        """
        pbsScript = dedent(pbsScript).format(name=directory,
                maxIdx=numJobs - 1, command = " ".join(cmd))
        subprocess.Popen(["qsub", "-"], stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
        qsub_stdout, qsub_stderr = subprocess.communicate(pbsScript)
        if qsub_stderr != "":
            print "WARNING: nonempty stderr response from qsub: ", qsub_stderr
        return qsub_stdout
