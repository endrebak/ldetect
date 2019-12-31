import pandas as pd
df = pd.read_table("/mnt/work/endrebak/ldetect/1kg/partition_covariances/CEU/chr22/16050075_16847854.tsv_old.gz", header=None, compression=None, sep="\s+", index_col=[2, 3])
df2 = pd.read_table("/mnt/work/endrebak/ldetect/1kg/partition_covariances/CEU/chr22/16050075_16847854.tsv.gz", header=None, compression=None, sep="\s+", index_col=[2, 3])

diff = df.index.difference(df2.index)
# print(df.head(10))
# r = df.iloc[:, 2:4].groupby(2).head(1)
# print(r[r[2] != r[3]])
