# Merges the second table into the first, printing the result to stdout.
# ***NOTE*** This means the initial stress is taken from the **first** table.
import pandas as pd
import sys
import csv

origTablePath, additionalTablePath = sys.argv[1:]

origTable       = pd.read_table(origTablePath)
additionalTable = pd.read_table(additionalTablePath)

origTable['file']       =       origTable['path'].replace(to_replace='.*/', regex=True, value=' ', inplace=False)
additionalTable['file'] = additionalTable['path'].replace(to_replace='.*/', regex=True, value=' ', inplace=False)

concat = origTable.append(additionalTable)
# always
grouped = concat.groupby("file", as_index=False)

merged = pd.DataFrame(columns=origTable.columns)
for name, group in grouped:
    newIdx = len(merged)
    if (len(group) == 1):
        merged.loc[newIdx] = group.iloc[0]
    elif (len(group) == 2):
        a = group.iloc[0]
        b = group.iloc[1]
        merged.loc[newIdx, :] = a;
        if (b['optimized'] < a['optimized']):
            merged.loc[newIdx, :] = b;
            if abs(a['targetE' ] - b['targetE' ] > 1e-8): raise Exception("Differing target E: " + str(a['targetE']) + ' and ' + str(b['targetE']))
            if abs(a['targetNu'] - b['targetNu'] > 1e-8): raise Exception("Differing target E: " + str(a['targetNu']) + ' and ' + str(b['targetNu']))
            # overwrite  all fields but initial stress level with b.
            # recompute reduction from this initial stress level
            merged.ix[newIdx, 'initial']   = a['initial'];
            reduction = b['optimized'] / a['initial']
            if (reduction > a['reduction']): raise Exception('Reduction factor worsened in merge')
            # if (reduction > b['reduction']): raise Exception('Reduction factor worsened in merge') # Note: this is actually allowed to happen because initial stress in a is treated as ground truth while b's init could be lower
            merged.ix[newIdx, 'reduction'] = b['optimized'] / a['initial']
    else: raise Exception("More than two rows grouped");

merged.drop('file', 1, inplace=True)
print merged.to_csv(sep="\t", index=False, quoting=csv.QUOTE_NONNUMERIC)
