import pandas as pd

partitions = snakemake.input["partitions"]
dfs = []
for p in partitions:
    df = pd.read_csv(p, sep=" ")
    dfs.append(df)

