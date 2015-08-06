from autocover_tools import autocoverRound, analyzeRuns
import json, os, sys, paths

if (len(sys.argv) != 3)
    print "usage: ./autocover_hpc.py roundNum config"
    print "example: ./autocover_hpc.py 1 2D_autocover_config.json"
    sys.exit(-1)

roundNum = int(sys.argv[1])
configPath = sys.argv[2]

config = json.load(file(configPath))

# If we've already run the full number of requested iterations, analyze the
# final iteration and quit
if (roundNum == maxIters + 1):
    analyzeRuns(None, roundNum, pat)
    return

# Otherwise, we launch this round's array job and then submit the next round as
# a dependent job.
arrayJobId = autocoverRound(roundNum, config['dim'], config['pattern'],
        paths.material(config['material']), config['targetERange'],
        config['targetNuRange'], config['targetNSubdiv'], hpc=True)

pbsScript = """\
###-----PBS Directives Start-----###

#PBS -V
#PBS -S /bin/bash
#PBS -N autocover_r{roundNum}
#PBS -l nodes=1:ppn=1
#PBS -l walltime=0:10:00
#PBS -l mem=4GB
#PBS -M fjp234@nyu.edu
#PBS -m a
#PBS -e localhost:${{PBS_O_WORKDIR}}/autocover_r{roundNum}.e${{PBS_JOBID}}
#PBS -o localhost:${{PBS_O_WORKDIR}}/autocover_r{roundNum}.o${{PBS_JOBID}}

###-----PBS Directives End-----###
cd ${{PBS_O_WORKDIR}}
{autocover} {roundNum} {config}
"""
pbsScript = dedent(pbsScript).format(
        autocover=os.path.dirname(os.path.realpath(__file__)),
        roundNum=roundNum + 1, config=configPath)

subprocess.Popen(["qsub", "-", "-Wafteranyarray:%s" % arrayJobId],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
qsub_stdout, qsub_stderr = subprocess.communicate(pbsScript)
print "qsub stdout: ", qsub_stdout
print "qsub stderr: ", qsub_stderr
