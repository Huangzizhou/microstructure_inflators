import sys, re
import pandas as pd

analyzedDataFile, = sys.argv[1:]

data = pd.read_table(analyzedDataFile)

for index, row in data.iterrows():
    pathStr = row['path']
    try: m = re.search('.*stdout_([0-9]+)_([0-9]+).txt', pathStr)
    except:
        print row
        print index
        raise Exception('Failed to apply regex to pathStr: ' + str(pathStr) + ' at index ' + str(index));
    if (not m): raise Exception('Invalid pathStr: ' + pathStr)
    pat,num = m.group(1), m.group(2)
    params=row['params']
    young=row['targetE']
    poisson=row['targetNu']
    print '{{"cwd": "/scratch/fjp234/wcsmin/autocover_coarser_p8_refit", "cmd": "/home/fjp234/microstructures/worst_case_stress/experiments/min_autocover/run_tensor_fit_from_params.sh {pat} {num} \'{params}\'", "stdout":"/scratch/fjp234/wcsmin/autocover_coarser_p8_refit/stdout_{pat}_{num}.txt", "stderr":"/scratch/fjp234/wcsmin/autocover_coarser_p8_refit/stderr_{pat}_{num}.txt"}},'.format(
            pat=pat, num=num, params=params)
