#!/usr/bin/env python
# Run a single autocoverage iteration. Sometimes this is a full "round," but
# because the round may require more jobs than can fit on the queue, often we
# must split rounds into sub-iterations:
#   The round "i" is started by passing i as the first argument. This creates
#   all the jobs files and launches up to config[maxSimultaneousJobs] of them.
#   If more remain, a autocover_hpc.py instance is submitted with
#   "i:maxSimultaneousJobs:numJobsInRound" as the first argument. This process
#   continues until all jobs for the round are submitted. Then an
#   autocover_hpc.py instance is submitted with "i+1" as the first argument.
from autocover_tools import *
from textwrap import dedent
from glob import glob
import json, os, sys, paths, subprocess

def usage():
    print "usage: ./autocover_hpc.py roundNum[:offset:numJobsInRound] config"
    print "example: ./autocover_hpc.py 1 2D_autocover_config.json"
    sys.exit(-1)

if len(sys.argv) != 3: usage()

# Decode roundNum[:subRoundJobOffset:numJobsInRound]
subRoundJobOffset, numJobsInRound = 0, 0
roundArg = sys.argv[1].split(':')
if (len(roundArg) not in [1, 3]): usage()
roundNum = int(roundArg[0])
if (len(roundArg) == 3):
    subRoundJobOffset, numJobsInRound = map(int, roundArg[1:])
    if (subRoundJobOffset <= 0): raise Exception("subRoundJobOffset must be strictly positive")
    
configPath = sys.argv[2]

config = json.load(file(configPath))
maxSimultaneousJobs = config.get('maxSimultaneousJobs', 500)

# If we've already run the full number of requested iterations, analyze the
# final iteration and quit
numIters = config['numIters']
if (roundNum == numIters + 1):
    outLUT = analyzeRuns(config['dim'], None, numIters, config['pattern'])
    outLUT.write(roundLUTPath(numIters))
    sys.exit(0)

# A subRoundJobOffset of 0 means this is the first time
# through--create the round's jobs and launch as many as possible.
opt = None
if (subRoundJobOffset == 0):
    opt = autocoverRoundOptimizer(roundNum, config)
    opt.writeJobs(roundDirectory(roundNum))
    numJobsInRound = len(opt.jobs)
else:
    # create dummy optimizer used to submit array jobs for subsequent sub-rounds
    opt = PatternOptimization(config)

if (subRoundJobOffset >= numJobsInRound): raise Exception("Error: invalid sub-round job offset")
numToSubmit = min(numJobsInRound - subRoundJobOffset, maxSimultaneousJobs)
nextOffset = subRoundJobOffset + numToSubmit
arrayJobId = opt.submitArrayJob(roundDirectory(roundNum), subRoundJobOffset, nextOffset - 1)

# Construct argument to pass the next autocover_hpc.py instance
nextRoundString = str(roundNum + 1)
if (nextOffset < numJobsInRound):
    nextRoundString = "%i:%i:%i" % (roundNum, nextOffset, numJobsInRound)

# NOTE: for some reason the -W depend... line cannot be passed to qsub
# via the subprocess arguments. It must be put in the batch script.
pbsScript = """\
###-----PBS Directives Start-----###

#PBS -V
#PBS -W depend=afteranyarray:{arrayJobId}
#PBS -S /bin/bash
#PBS -N {name}
#PBS -l nodes=1:ppn=1
#PBS -l walltime=0:10:00
#PBS -l mem=4GB
#PBS -M fjp234@nyu.edu
#PBS -m a
#PBS -e localhost:${{PBS_O_WORKDIR}}/{name}.e${{PBS_JOBID}}
#PBS -o localhost:${{PBS_O_WORKDIR}}/{name}.o${{PBS_JOBID}}

###-----PBS Directives End-----###
cd ${{PBS_O_WORKDIR}}
{autocover} {roundString} {config}
"""
pbsScript = pbsScript.format(name="ac_%i_r%s" % (config['pattern'], nextRoundString),
        autocover=os.path.realpath(__file__),
        pat=config['pattern'],
        roundString=nextRoundString, config=os.path.realpath(configPath),
        arrayJobId=arrayJobId)
p = subprocess.Popen(["qsub", "-"],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
qsub_stdout, qsub_stderr = p.communicate(pbsScript)
print "qsub stdout: ", qsub_stdout,
print "qsub stderr: ", qsub_stderr
