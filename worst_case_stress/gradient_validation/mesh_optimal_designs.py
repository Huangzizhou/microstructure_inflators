import pandas as pd
import subprocess
st2 = pd.HDFStore('d2_stats.hdf')
for p in ['inf', '12.0']:
    t = st2['minimizer']['0.00004'].transpose()[p]
    for k in t.keys():
        subprocess.call(['../LpHoleInflator_cli', str(k), str(t[k]), 'L%s_r%s.msh' % (p, k), '-a', '0.00004'])
st2.close()
