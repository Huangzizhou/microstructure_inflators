import subprocess, numpy as np
import os
for i, b in enumerate(np.arange(1.0, 0.25, -0.01)):
    subprocess.call(['./ellipse_newton', '1', str(b)])
    subprocess.call(['gnuplot', 'plot.gpi'])
    os.rename('iterations.png', 'iterations_%04d.png' % i)
    # os.rename('backtracking_iterations.png', 'backtracking_%04d.png' % i)
    # os.rename('angle.png', 'angle_%04d.png' % i)
    os.rename('cp_x.png', 'cp_x_%04d.png' % i)
    os.rename('distance.png', 'distance_%04d.png' % i)
subprocess.call(['images_to_mp4.sh', '.', 'iterations.mp4', 'iterations_%04d.png'])
subprocess.call(['images_to_mp4.sh', '.', 'cp_x.mp4', 'cp_x_%04d.png'])
# subprocess.call(['images_to_mp4.sh', '.', 'backtracking.mp4', 'backtracking_%04d.png'])
# subprocess.call(['images_to_mp4.sh', '.', 'angle.mp4', 'angle_%04d.png'])
subprocess.call(['images_to_mp4.sh', '.', 'distance.mp4', 'distance_%04d.png'])
