import pandas as pd
import numpy as np


def update_covar_and_loci(covar, loci, pos_to_delete_below, partition_file):

    covar = covar[covar.i >= pos_to_delete_below]
    loci = loci[loci >= pos_to_delete_below]

    to_add = read_partition(partition_file)

    covar = covar.append(to_add).reset_index(drop=True)
    loci = loci.append(to_add.i.drop_duplicates()).reset_index(drop=True)

    g = covar.i.reset_index().groupby("i").index
    covar_starts = g.first().values
    covar_ends = g.last().values + 1

    return covar, loci, covar_starts, covar_ends


def find_start_locus(curr_locus, loci, cutoff):

    if curr_locus < 0:
        for i, locus in enumerate(loci):
            if locus >= cutoff:
                return i, locus
    else:
        return loci.index(curr_locus), curr_locus


def find_end_locus(p_num, partitions, loci, cutoff):

    if p_num + 1 < len(partitions):
        end_locus = partitions.iloc[p_num + 1, 0] 
        end_locus_index = -1
        return end_locus_index, end_locus
    else:
        for i in reversed(range(0, len(loci))):
            if loci[i] <= cutoff:
                return i, loci[i]


def read_partitions(pfiles):

    if pfiles:
        df = pd.concat([read_partition(f) for f in pfiles]).reset_index(drop=True)
        return df, df.i.drop_duplicates().reset_index(drop=True) 
    else:
        df = pd.DataFrame({"i": np.array([], dtype=int), "j": np.array([], dtype=int), "val": np.array([], dtype=float)})
        return df, df.i


def read_partition(pfile):

    if pfile.endswith(".pq"):
        return pd.read_parquet(pfile)
    else:
        return pd.read_table(pfile, sep=" ", usecols=[2, 3, 7], names="i j val".split())


def covar_files_map(covariance_files, partitions):

    """Map of interval to covariance file."""

    covar_dict = {}
    for f, (start, end) in zip(covariance_files, partitions.itertuples(index=False)):
        covar_dict[start, end] = f

    return covar_dict


def files_to_preread(partitions, covar_dict, larger_than):

    i = 0
    to_preread = []
    while i < len(partitions) - 1 and larger_than >= partitions[i + 1][0]:
        pstart, pend = partitions[i]
        to_preread.append(covar_dict[pstart, pend])
        i += 1

    return to_preread, i


def preread_files(partitions, covar_dict, larger_than):

    to_preread, i = files_to_preread(partitions, covar_dict, larger_than)
    covar, loci = read_partitions(to_preread)

    return covar, loci, i
