import pandas as pd
import numpy as np
from ldetect2.src.m2v import mat2vec

from sys import argv

partitions = argv[1]
theta2 = argv[2]

covariances = argv[3:]

theta2 = float(open(theta2).readline().strip())

import sys

# print(partitions, file=sys.stderr)
# print(covariances, file=sys.stderr)

# partitions = snakemake.input["partitions"]
# covariances = snakemake.input["covariances"]




dfs = []
from time import time
length = 0
memory_use = 0
for i, c in enumerate(covariances):
    start = time()
    df = pd.read_parquet(c)
    length += len(df)
    memory_use += df.memory_usage(deep=True).sum()
    dfs.append(df)
    end = time()
    print(i, c, i/len(covariances), end - start, length, memory_use / 1e9, file=sys.stderr)


covariances = sorted(covariances, key=lambda k: k.split("/"))

df = pd.concat(dfs)

print("Done concatenating!")
print("df.memory_usage(deep=True).sum()", df.memory_usage(deep=True).sum())

# s = pd.read_csv(covariances[-1], sep=" ", usecols=[2], names=["i"],
#                  dtype={"i": np.int32}, squeeze=True)
# max_ = s.max()

ps = pd.read_table(partitions, sep=" ", header=None)
new_ends = ((ps[0] + ps[1].shift(-1).values)/2)
new_ends = new_ends.fillna(df.i.max()).astype(np.int32)
ps.insert(ps.shape[1], 2, new_ends)

# assert len(ps) == len(covariances), "Number of partitions and covariance files are not the same!"

mat2vec(df.i.values, df.j.values, df.val.values, ps, theta2)
