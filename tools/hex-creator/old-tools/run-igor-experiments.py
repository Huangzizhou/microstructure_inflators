#!/usr/bin/env python
import os
from subprocess import call

cwd = os.getcwd()

if os.path.isfile('igor-instances/negative-poisson_p1-0.97_p2-60_p3-0.94_p4-0.9.msh'):
    print "Already computed!"
else:
    cmd = [cwd + '/auxetic-chiral-creator.py', '0.97', '60', '0.94', '0.9', 'igor-instances/negative-poisson_p1-0.97_p2-60_p3-0.94_p4-0.9.wire', 'igor-instances/negative-poisson_p1-0.97_p2-60_p3-0.94_p4-0.9.msh']
    call(cmd)

if os.path.isfile('igor-instances/negative-poisson_p1-0.975_p2-30_p3-0.935_p4-0.72.msh'):
    print "Already computed!"
else:
    cmd = [cwd + '/auxetic-chiral-creator.py', '0.975', '30', '0.935', '0.72', 'igor-instances/negative-poisson_p1-0.975_p2-30_p3-0.935_p4-0.72.wire', 'igor-instances/negative-poisson_p1-0.975_p2-30_p3-0.935_p4-0.72.msh']
    call(cmd)

if os.path.isfile('igor-instances/negative-poisson_p1-0.97_p2-30_p3-0.925_p4-0.86.msh'):
    print "Already computed!"
else:
    cmd = [cwd + '/auxetic-chiral-creator.py', '0.97', '30', '0.925', '0.86', 'igor-instances/negative-poisson_p1-0.97_p2-30_p3-0.925_p4-0.86.wire', 'igor-instances/negative-poisson_p1-0.97_p2-30_p3-0.925_p4-0.86.msh']
    call(cmd)

if os.path.isfile('igor-instances/negative-poisson_p1-0.99_p2-30_p3-0.94_p4-0.6.msh'):
    print "Already computed!"
else:
    cmd = [cwd + '//auxetic-chiral-creator.py', '0.99', '30', '0.94', '0.6', 'igor-instances/negative-poisson_p1-0.99_p2-30_p3-0.94_p4-0.6.wire', 'igor-instances/negative-poisson_p1-0.99_p2-30_p3-0.94_p4-0.6.msh']
    call(cmd)


if os.path.isfile('igor-instances/positive-poisson_p1-0.9_p2-25_p3-1.0_p4-0.6.msh'):
    print "Already computed!"
else:
    cmd = [cwd + '/hexa-many-pillars-igor-parameters.py', '0.9', '25', '1.0', '0.6', 'igor-instances/positive-poisson_p1-0.9_p2-25_p3-1.0_p4-0.6.wire', 'igor-instances/positive-poisson_p1-0.9_p2-25_p3-1.0_p4-0.6.msh']
    call(cmd)

if os.path.isfile('igor-instances/positive-poisson_p1-0.87_p2-25_p3-1.0_p4-0.6.msh'):
    print "Already computed!"
else:
    cmd = [cwd + '/hexa-many-pillars-igor-parameters.py', '0.87', '25', '1.0', '0.6', 'igor-instances/positive-poisson_p1-0.87_p2-25_p3-1.0_p4-0.6.wire', 'igor-instances/positive-poisson_p1-0.87_p2-25_p3-1.0_p4-0.6.msh']
    call(cmd)

if os.path.isfile('igor-instances/positive-poisson_p1-0.87_p2-25_p3-1.0_p4-0.8.msh'):
    print "Already computed!"
else:
    cmd = [cwd + '/hexa-many-pillars-igor-parameters.py', '0.87', '25', '1.0', '0.8', 'igor-instances/positive-poisson_p1-0.87_p2-25_p3-1.0_p4-0.8.wire', 'igor-instances/positive-poisson_p1-0.87_p2-25_p3-1.0_p4-0.8.msh']
    call(cmd)

if os.path.isfile('igor-instances/positive-poisson_p1-0.93_p2-25_p3-1.0_p4-0.7.msh'):
    print "Already computed!"
else:
    cmd = [cwd + '/hexa-many-pillars-igor-parameters.py', '0.93', '25', '1.0', '0.7', 'igor-instances/positive-poisson_p1-0.93_p2-25_p3-1.0_p4-0.7.wire', 'igor-instances/positive-poisson_p1-0.93_p2-25_p3-1.0_p4-0.7.msh']
    call(cmd)