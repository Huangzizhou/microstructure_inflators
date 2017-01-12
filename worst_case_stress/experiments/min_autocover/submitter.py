#!/usr/bin/env python
# Splits up an array job into bits that can be placed on the queue.
from getpass import getuser
import sys, re, subprocess, time
user = getuser()

MAX_JOBS = 350;
MAX_PER_SUBMIT = 200;

# Path to the .pbs script
arrayJobPBS, = sys.argv[1:]

arraySize = None
for l in file(arrayJobPBS):
    m = re.match('#PBS -t 0-([0-9]+)', l)
    if m:
        if (arraySize != None): raise Exception("Multiple array size lines encountered");
        arraySize = int(m.group(1))

numSubmitted = 0;

def submitMore():
    global numSubmitted # modified in this function

    qstatOut = subprocess.check_output(['qstat', '-t', '-u', user])
    held = []
    running = []
    queued = []
    for l in qstatOut.strip().split('\n'):
        m = re.match('([0-9\[\]]+)\s+.*\s([HRQ])\s', l)
        if m:
            jobid = m.group(1)
            stat = m.group(2)
            if   (stat == 'H'):    held.append(jobid)
            elif (stat == 'R'): running.append(jobid)
            elif (stat == 'Q'):  queued.append(jobid)
            else: raise Exception("Unknown status " + stat)
    numInQueue = len(held) + len(running) + len(queued)
    numSubmitNow = min(min(arraySize - numSubmitted, MAX_PER_SUBMIT),
                       MAX_JOBS - numInQueue)

    if (numSubmitNow <= 0): return

    try:
        subprocess.check_call(['qsub', arrayJobPBS, '-t',
                               str(numSubmitted) + '-' +
                               str(numSubmitted + numSubmitNow - 1)])
        numSubmitted += numSubmitNow
    except:
        print "WARNING: Bad return value from qsub; assuming jobs not submitted";

while (numSubmitted < arraySize):
    try: submitMore()
    except Exception as e: print "EXCEPTION: ", e
    time.sleep(600) # every 10 minutes
