#!/usr/bin/env python
# Monitor this user's current jobs for held jobs that are "stuck."
# Torque has a bug where dependent jobs can be held from the queue
# forever, even after their dependencies have finished. We detect this and
# release the job to the queue.
#
# Also, for intermediate sub-rounds (which only launch new jobs and don't
# actually depend on the previous jobs' results), we release dependent jobs
# into the queue when all of their dependencies are already running and fewer
# than <DEP_MIN> remain. This avoids the situation where a few slow jobs hold up
# future batches. There's a small chance of starting the next round (analyzing
# job output) before all jobs from previous subrounds have finished,
# but that's fine--this next round will try for the missing gridpoints again.
import subprocess, re, time
from getpass import getuser
user = getuser()

DEP_MIN=1

def releaseStuckDependents():
    qstatOut = subprocess.check_output(['qstat', '-t', '-u', user])
    held = []
    running = []
    queued = []
    for l in qstatOut.strip().split('\n'):
        m = re.match('([0-9\[\]]+)\s+.*\s([HRQ])\s', l)
        if m:
            jobid = m.group(1)
            stat = m.group(2)
            if (stat == 'H'):      held.append(jobid)
            elif (stat == 'R'): running.append(jobid)
            elif (stat == 'Q'):  queued.append(jobid)
            else: raise Exception("Unknown status " + stat)

    # Get the ID part of a job only (for array jobs, this is a substring up to '[')
    def baseJobId(j): return j.split('[')[0]

    # print "currently held: ", "\t".join(held)
    for h in held:
        fullName = subprocess.check_output('qstat -f %s | grep Job_Name' % h, shell=True).split('=')[1].strip()
        depend = ""
        try: depend = subprocess.check_output('qstat -f %s | grep depend' % h, shell=True).strip()
        except: pass
        if (depend == ""): raise Exception("Empty dependency for job %s (qstat glitch--job may eventually run)" % h)
        # Currently we only support a single (possibly array) dependency
        m = re.match('depend = (after[^:]*):([0-9]+)\[?\]?.soho\S*$', depend)
        if (not m): raise Exception("Unrecognized dependency line")
        depType = m.group(1)
        depJobId = m.group(2)

        queuedDepCount, runningDepCount = 0, 0
        for qj in queued:
            if baseJobId(qj) == depJobId: queuedDepCount += 1
        for rj in running:
            if baseJobId(rj) == depJobId: runningDepCount += 1
        intermediateJob = len(fullName.split(':')) == 3
        if (intermediateJob and (queuedDepCount == 0) and (runningDepCount < DEP_MIN)):
            print "%s (%s) depends on %i queued and %i running jobs... " % (h, fullName, queuedDepCount, runningDepCount)
            qaCommand = ['qalter', '-W', 'depend=%s' % depType, h]
            print " ".join(qaCommand)
            subprocess.call(qaCommand)

while True:
    try: releaseStuckDependents()
    except Exception as e: print "EXCEPTION: ", e
    time.sleep(300)
