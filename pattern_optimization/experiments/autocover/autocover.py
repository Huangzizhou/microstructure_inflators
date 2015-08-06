import json, sys, paths
from autocover_tools import autocoverRound

if (len(sys.argv) != 3):
    print "usage: ./autocover.py startingRoundNumber config"
    print "example: ./autocover.py 1 2D_autocover_config.json"
    sys.exit(-1)

roundNum = int(sys.argv[1])
configPath = sys.argv[2]

config = json.load(file(configPath))

numIters = config['numIters']
for roundNum in range(roundNum, numIters):
    print "Running round %i/%i" % (roundNum, numIters)
    autocoverRound(roundNum, config['dim'], config['pattern'],
            paths.material(config['material']), config['targetERange'],
            config['targetNuRange'], config['targetNSubdiv'])
