from glob import glob
files = glob("/mnt/work/endrebak/ldetect/1kg/partition_covariances/CEU/chr22/*")

import pandas as pd
import numpy as np

d = {}

dfs = []
for i, f in enumerate(files):
    print(i/len(files))
    df = pd.read_csv(f, sep=" ", usecols=[2, 3, 7], names="i j val".split(), dtype={"i": np.int32, "j": np.int32})
    dfs.append(df)
    # for k1, k2, v in df.itertuples(index=False):
    #     d[k1, k2] = v

df = pd.concat(dfs)
# print(len(d))


