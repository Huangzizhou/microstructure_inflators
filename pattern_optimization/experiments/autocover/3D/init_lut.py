import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
import LookupTable

# Usage: python init_lut.py pattern out_path
pattern,path = sys.argv[1:]
pattern = int(pattern)
l = LookupTable.extract(3, pattern, '.', init=True) 
l.write(path)
