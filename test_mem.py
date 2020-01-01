from glob import glob
files = glob("/mnt/work/endrebak/ldetect/1kg/partition_covariances/CEU/chr22/*.pq")

import pandas as pd
import numpy as np

d = {}

dfs = []
for i, f in enumerate(files):
    print(i/len(files))
    df = pd.read_parquet(f)
    dfs.append(df)
    # for k1, k2, v in df.itertuples(index=False):
    #     d[k1, k2] = v

# df = pd.concat(dfs)
# print(len(d))


