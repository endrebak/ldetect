import pandas as pd
import numpy as np
from ldetect2.src.matrix_to_vector import mat2vec

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


ps = pd.read_table(partitions, sep=" ", header=None)

dfs = []
for c in covariances:
    # print(c)
    df = pd.read_csv(c, sep=" ", usecols=[2, 3, 7], names="i j val".split())
    dfs.append(df)

df = pd.concat(dfs)

# print("before")
mat2vec(df.i.astype(np.int32).values, df.j.astype(np.int32).values, df.val.values, ps, theta2)
# print("after")
