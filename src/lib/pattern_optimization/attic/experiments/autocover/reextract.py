# Re-extract the round lookup tables from archived stdout files.
# (Useful for changing the filtering thresholds/parameters).
import sys, os
import os.path
from glob import glob
import tarfile, shutil
import numpy as np
from LookupTable import LUT, extract as extractLUT
dim, pat, rounds_dir, = sys.argv[1:]

dim = int(dim)
pat = int(pat)

os.chdir(rounds_dir)
try: shutil.rmtree('reextract')
except: pass

os.mkdir('reextract')

# np.seterr(all='raise')

prevLUT = LUT("round_0000.txt")
for xiv in sorted(glob('*.tgz')): # visit rounds in order so that we can concatenate
    archive = tarfile.open(xiv, 'r')
    names = archive.getnames()
    extractDir = None
    for n in sorted(names):
        if n[-4:] == '.txt':
            if extractDir == None: extractDir = os.path.dirname(n)
            elif extractDir != os.path.dirname(n): raise Exception('Mismatch directory names')
            f = archive.extract(n, 'reextract')

    lut = extractLUT(dim, pat, 'reextract/' + extractDir, printableOnly=True)
    lut.filterAnisotropy(0.99, 1.01)

    lut.union(prevLUT)
    prevLUT = lut

    lut.write('reextract/' + extractDir + '.txt')
    try: shutil.rmtree('reextract/' + extractDir)
    except: pass
