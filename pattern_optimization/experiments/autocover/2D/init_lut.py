import sys
sys.path.append('..')
import LookupTable
# Usage: python init_lut.py pattern out_path
pattern,path = sys.argv[1:]
pattern = int(pattern)
l = LookupTable.extract(2, pattern, '.', init=True) 
l.write(path)
