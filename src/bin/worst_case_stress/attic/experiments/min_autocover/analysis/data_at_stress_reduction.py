import sys, pandas as pd, numpy as np

dataTablePath,reduction = sys.argv[1:]

data = pd.read_table(dataTablePath)

# Data table reduction field is opt/initial, whereas the command line argument
# is initial/opt.
indices = np.where(data['reduction'] < 1.0 / float(reduction) + 1e-9)[0]

print data.ix[indices].to_csv(sep="\t", index=False)
