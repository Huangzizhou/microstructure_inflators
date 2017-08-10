import json

jobs = json.load(open("jobs.json"))
for j in jobs:
    sweepType=j["cwd"].split('/')[-1]
    print "{} > {}".format(j["cmd"], sweepType + "_" + j["stdout"])
