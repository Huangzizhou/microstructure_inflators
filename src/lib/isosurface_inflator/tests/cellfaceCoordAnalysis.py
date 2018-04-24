import sys, numpy as np

coordPath,coordIdx,faceCoord = sys.argv[1:]
coordIdx = int(coordIdx)
faceCoord = float(faceCoord)

coords = [map(float, l.strip().split()) for l in file(coordPath)]

result = [abs(c[coordIdx] - faceCoord) for c in coords]
sortIdx = np.argsort(result)
for i in sortIdx:
    print "{}: {}".format(coords[i], result[i])
