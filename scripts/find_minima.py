import pandas as pd

from sys import argv
import scipy.ndimage.filters as filters
import scipy.signal as sig

from sklearn.utils.fast_dict import IntFloatDict as fdict

def create_fastdict(keys=None, vals=None):
    if keys is None:
        d = fdict(np.array([], dtype=int), np.array([]))
        return d
    else:
        return fdict(keys, vals)

import math
sqrt = math.sqrt
import numpy as np

from time import time, sleep

f = argv[1]

df = pd.read_table(f, header=None, names="pos val".split())

partition_file = argv[2]
# print(partition_file)
partitions = pd.read_table(partition_file, sep=" ", header=None)
partitions.index = range(len(partitions))
partitions.columns = ["Start", "End"]
# print(partitions)

covariance_files = argv[3:]


def _apply_filter(arr, width):

    "should probably memoize this function :"

    val = (2 * width + 1)
    moving_avg_a = np.ones(val) * 1/(val)
    # start = time()
    # print("get_window")
    a = sig.get_window('hanning', val)

    # print("convolve1d")
    ga = filters.convolve1d(arr, a/a.sum())

    # print("argrelextrema")
    minima_a = sig.argrelextrema(ga, np.less)[0]

    return minima_a

def _find_end(a, init_search_location, n_bpoints):
    search_location = init_search_location
    # print("a", a)
    filtered = _apply_filter(a, search_location)
    # print("filtered", filtered)
    while len(filtered) >= n_bpoints:
        # print("search_location", search_location)
        filtered = _apply_filter(a, search_location * 2)
        search_location = search_location * 2

    return search_location

# "When looking up a index, you use that index as the width for the filter"

# "binary search but apply the filter for each time"

def _width_raw(a, n_bpoints):

    "gives the same result as the variable found_width in the original code"

    search_val = n_bpoints

    filtered = _apply_filter(a, search_val)
    idx = np.searchsorted(filtered, search_val, side="left") - 1

    return len(filtered) - idx


def _trackback(array, search_val, start_search):

    step_coarse = 20
    delta_coarse = 200
    step_fine = 1

    for i in range(start_search + step_coarse, len(array), step_coarse):

        filtered = _apply_filter(array, i)

        if len(filtered) == search_val:
            start_search = i


    for i in range(start_search + step_fine, len(array), step_fine):

        filtered = _apply_filter(array, i)

        if len(filtered) == search_val:
            start_search = i

    return start_search

def _width(val, n_bpoints):

    search_value = 1000
    end = _find_end(val, search_value, n_bpoints)

    found_width = end - _width_raw(val, n_bpoints)

    found_width_trackback_raw = end - _trackback(val, search_value, found_width)

    return found_width_trackback_raw


def find_breakpoint_loci(df, n_snps=50):

    n_bpoints = int(math.ceil(len(df) / n_snps) - 1)

    val = df.val
    pos = df.pos.values

    width = _width(val, n_bpoints)

    filtered = _apply_filter(val, width)
    return pos[filtered]


def compute_zero_metric(loci, partitions, breakpoints):

    # end locus is just start of next
    # or last snp smaller than last partition end

    # iterate over the partitions
    # for each locus, check if the current locus is beyond the current breakpoint
    #   if so, calculate the statistics for the breakpoint

    # afterwards for each of the snps correlated with curr_locus
    # if they are above the breakpoint index, add the correlation coefficient (x, y / ((x, x) * (y, y)))

    # for each iteration, increase block_width with one
    #                     total_N_SNPs with one

    # can compute N_zero and N_nonzero on first pass
    #   also, store all corr_coeffs that should be computed

    curr_locus_index = 0
    block_height = 0
    block_width = 0
    block_width_sum = 0
    total_snps = 0
    nzero = 0
    # breakpoints_iter = iter(breakpoints)
    curr_breakpoint_index = 0
    curr_breakpoint = breakpoints[curr_breakpoint_index]

    # with open("curr_locus_and_breakpoints.txt", "w+") as clb:
    #     for start, end in partitions.itertuples(index=False):
    loci_to_compute_later = np.zeros(len(loci), dtype=int)
    breakpoints_out = np.zeros(len(loci), dtype=int)


    for i in range(len(loci)):
        curr_locus = loci[curr_locus_index]
        # while curr_locus < end:
        curr_breakpoint = breakpoints[curr_breakpoint_index]
        if curr_locus > curr_breakpoint:

            block_height = 0 - total_snps
            nzero += block_height * block_width
            block_width_sum += block_width

            # print("-----" * 5)
            # print("curr_locus", curr_locus)
            # print("curr_breakpoint", curr_breakpoint)
            # print("block_height", block_height)
            # print("block_width", block_width)
            # print("nzero", nzero)
            # print("block_width_sum", block_width_sum)

            block_width = 0

            curr_breakpoint_index += 1

        if curr_breakpoint_index >= len(breakpoints):
            break

        # print(total_snps)
        loci_to_compute_later[total_snps] = curr_locus
        breakpoints_out[total_snps] = breakpoints[curr_breakpoint_index]

        block_width += 1

        if curr_locus_index + 1 < len(loci):
            curr_locus_index += 1
            curr_locus = loci[curr_locus_index]
            total_snps += 1

    nzero += total_snps * block_width_sum

    # print("len(loci)", len(loci))
    return nzero, loci_to_compute_later[:total_snps], breakpoints_out[:total_snps]


def compute_sum_and_nonzero(loci, bps, i_, j_, covars, autocovar):

    i = 0
    j = 0
    nonzero = 0
    metric_sum = 0
    covar_len = len(covars)
    for i in range(len(loci)):
        locus = loci[i]
        breakpoint = bps[i]

        # print("i", i)
        # print(covars.iloc[j])

        while j < covar_len and i_[j] != locus:
            # print("!= locus", j)
            j += 1

        while j < covar_len and i_[j] == locus:
            # print("== locus", j)
            if j_[j] > breakpoint:
                corrcoeff = covars[j] / sqrt(autocovar[i_[j]] * autocovar[j_[j]])
                metric_sum += corrcoeff ** 2
                # print(i_[j], j_[j], covars[j], autocovar[i_[j]], autocovar[j_[j]], corrcoeff ** 2)
                nonzero += 1

            j += 1

             
        # print("j", j)

    return metric_sum, nonzero





def initialize_search(covars, autocovars):

    """TODO: ensure that covars only between start and ends of search

    For first block, include all snps, even though they are not above search start"""

    # from sklearn.utils.fast_dict import IntFloatDict as fdict

    iloci = covars.i.drop_duplicates()
    jloci = covars.j.drop_duplicates()
    length = len(iloci)
    # v, h = np.zeros(length), np.zeros(length)
    v, h = dict(zip(iloci, [0] * len(iloci))), dict(zip(jloci, [0] * len(jloci)))


    for i, j, val in covars.itertuples(index=False):
        corr_coeff = (val / math.sqrt(autocovars[i] * autocovars[j])) ** 2
        # print("ijval", i, j, val, corr_coeff, autocovars[i], autocovars[j])
        # if i == 39967768:
        #     print(i, j)
        v[i] += corr_coeff
        h[j] += corr_coeff

    # print(v[39967768])
    # sleep(2)
    # print(h[39967768])
    # raise

    return v, h
    # only compute between the breakpoints

    # while curr_locus <= end_locus:

    # for first block (initial_breakpoint) always include snps, otherwise
    # require that curr_locus > snp_first




    # horizontal = fdict()


# def local_search_right(v, h, partitions, metric_sum, zero, nonzero):





def search_starts_ends(breakpoints, partitions):


    diff = pd.Series(breakpoints).diff()
    # print(diff)
    search_starts = (breakpoints + (diff/2).shift(-1)).shift()
    search_starts[0] = partitions.iloc[0, 0]
    search_starts = search_starts.astype(int)
    search_ends = search_starts.shift(-1)
    search_ends.values[-1] = partitions.iloc[-1, 1]

    search_ends = search_ends.astype(int)

    bps = pd.Series(breakpoints)
    bps_next = bps.shift(-1)
    bps_next.values[-1] = partitions.iloc[-1, 1]
    bps_next = bps_next.astype(int)

    bps_last = bps.shift(1)
    bps_last.values[0] = partitions.iloc[0, 0]
    bps_last = bps_last.astype(int)
    df = pd.concat([search_starts, search_ends, bps, bps_last, bps_next], axis=1)
    df.columns = ["SearchStart", "SearchEnd", "Breakpoint", "LastBreakpoint", "NextBreakpoint"]

    # partition_end = relevant_partitions(partitions, )
    # search_end_last = relevant_partitions(partitions, )

    return df


def local_search_right(loci, v, h, breakpoint, snp_bottom_idx, snp_top_idx, metric_sum, zero_metric):

    idx = np.searchsorted(loci, breakpoint) + 1
    print("len(loci)", len(loci))
    # print("idx", idx, len(loci), print(type(loci)))
    # print("idx", idx, loci.values[idx])
    init_breakpoint_locus = loci[idx]
    loci = loci[idx:]
    # print("loci")
    # print(loci)

    n_horiz = 0
    n_vert = 0
    n_curr = zero_metric
    curr_sum = metric_sum
    min_metric = curr_sum / zero_metric
    min_breakpoint = -1
    min_distance_right = -1

    print("curr_loc", loci[0])
    for curr_loc in loci:

        curr_sum = curr_sum - h[curr_loc] + v[curr_loc]
        # print("curr_sum", curr_sum)
        # print("idx", idx)
        n_horiz = idx - snp_bottom_idx - 1
        # print("n_horiz", n_horiz)
        n_vert = snp_top_idx - idx
        # print("n_vert", n_vert)
        n_curr = n_curr - n_horiz + n_vert
        # print("n_curr", n_curr)

        curr_metric = curr_sum / n_curr
        # print(curr_loc, "curr_metric", curr_metric, h[curr_loc], v[curr_loc])
        # print("curr_loc", curr_loc, "curr_metric", curr_metric)
        # print("idx", idx)
        if curr_metric < min_metric:
            min_metric = curr_metric
            min_breakpoint = curr_loc
            min_distance_right = curr_loc - init_breakpoint_locus

        idx += 1

    # print("min_metric", min_metric, min_breakpoint, min_distance_right)

    return min_metric, min_breakpoint, min_distance_right

    # # find the index in locus_list less or equal to the current breakpoint
    # # that is the start of the search?
    # # raise

    

    # pass


def locs_to_use(covars, count, end, start, next_breakpoint):

    if count != 0:
        good_covars_idx = (covars.i < end) & (covars.j <= next_breakpoint) & (covars.i > start)
        good_covars = covars[good_covars_idx]
        good_loci_idx = (loci < end) & (loci > start)
        good_loci = loci[good_loci_idx]

    else:
        good_covars_idx = (covars.i < end) & (covars.j <= next_breakpoint) 
        good_loci = loci[loci < end]
        good_covars = covars[good_covars_idx]

    return good_covars, good_loci.values


    """
self.snp_first = start_search
self.snp_last = stop_search

tmp_partitions = flat.get_final_partitions(self.input_config, self.name, start_search, stop_search)

if initial_breakpoint_index+1 < len(breakpoints):
    self.snp_top = breakpoints[initial_breakpoint_index+1]
else:
    self.snp_top = tmp_partitions[len(tmp_partitions)-1][1]

# This is the bottom bound for the search space (bottom border)
if initial_breakpoint_index-1 >= 0:
    self.snp_bottom = breakpoints[initial_breakpoint_index-1]
else:
    self.snp_bottom = tmp_partitions[0][0]
    """

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


def find_start_locus(curr_locus, loci, snp_bottom):

    if curr_locus < 0:
        for i, locus in enumerate(loci):
            if locus >= snp_bottom:
                return i, locus
    else:
        return loci.index(curr_locus), curr_locus

def find_end_locus(p_num, partitions, loci, snp_last):

    if p_num + 1 < len(partitions):
        end_locus = partitions.iloc[p_num + 1, 0] # start of next partition
        end_locus_index = -1
        return end_locus_index, end_locus
    else:
        for i in reversed(range(0, len(loci))):
            if loci[i] <= snp_last:
                return i, loci[i]


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


# def precompute_values(covar_dict, snp_top, snp_bottom, snp_first, snp_last):

#     v = create_fastdict()
#     h = create_fastdict()

#     snp_top = 0
#     snp_bottom = 0
#     snp_first = 0
#     snp_last = 0
#     initial_breakpoint_index = 0
#     end_locus = 0

#     covar_dict_items = covar_dict.items()

#     for p_num, ((partition_start, partition_end), f) in enumerate(covar_dict_items):

#         to_preread = []
#         last_pnum = -1
#         if snp_bottom >= partition_start:
#             to_preread.append(f)
#             last_pnum = p_num
#         else:
#             break

#         covar, loci = read_partitions(to_preread)

#         curr_locus = -1
#         for p_num, ((partition_start, partition_end), f) in enumerate(covar_dict_items[p_num + 1], p_num + 1):

#             covar, loci, covar_starts, covar_ends = update_covar_and_loci(covar, loci, end_locus, f)
#             autocovars = covar[covar.i == covar.j]
#             autocovar = create_fastdict(autocovars.i.values, autocovars.val.values)

#             start_locus_index, start_locus = find_start_locus(curr_locus, loci, snp_bottom)
#             curr_locus_index, curr_locus = start_locus_index, start_locus
#             end_locus_index, end_locus = find_end_locus(p_num, curr_locus, loci, snp_last)

#             precomputed_loci = []

#             # for curr_locus_index in range(curr_locus_index, len(loci)):
#             while loci[curr_locus_index] <= end_locus:

#                 precomputed_loci.append(curr_locus)

#                 if (curr_locus > snp_first or initial_breakpoint_index == 0) and (curr_locus <= snp_last):

#                     covar_start, covar_end = covar_starts[curr_locus_index], covar_ends[curr_locus_index]

#                     for x in range(covar_start, covar_end):

#                         if covar[x][1] <= snp_top:

#                             i = covar[x][0]
#                             j = covar[x][1]
#                             val = covar[x][2]
#                             corr_coeff = (val / math.sqrt((autocovar[i] * autocovar[j]))) ** 2
#                             v[i] += corr_coeff
#                             h[j] += corr_coeff

#                 curr_locus_index += 1

#         return v, h


"""
self.snp_first = start_search
self.snp_last = stop_search

tmp_partitions = flat.get_final_partitions(self.input_config, self.name, start_search, stop_search)

if initial_breakpoint_index+1 < len(breakpoints):
    self.snp_top = breakpoints[initial_breakpoint_index+1]
else:
    self.snp_top = tmp_partitions[len(tmp_partitions)-1][1]

# This is the bottom bound for the search space (bottom border)
if initial_breakpoint_index-1 >= 0:
    self.snp_bottom = breakpoints[initial_breakpoint_index-1]
else:
    self.snp_bottom = tmp_partitions[0][0]
"""



def relevant_partitions(partitions, start, end):

    for i, (pstart, pend) in enumerate(partitions.itertuples(index=False)):

        if pstart <= start and pend >= start:
            xstart = i

        if pstart <= end and pend >= end:
            xend = i + 1

    rp = partitions.iloc[xstart:xend]

    return list(zip(rp.Start, rp.End))


def precompute_values(partitions, search_start, search_end, last_breakpoint, next_breakpoint):

    iter_partitions = relevant_partitions(partitions, search_start, search_end)

    curr_locus = -1

    i = 0
    to_preread = []
    while i < len(iter_partitions) - 1 and last_breakpoint >= iter_partitions[i][0]:
        pstart, pend = iter_partitions[i]
        to_preread.append(covar_dict[pstart, pend])
        i += 1

    covar, loci = read_partitions(to_preread)

    snp_bottom, snp_top, snp_first, snp_last = last_breakpoint, next_breakpoint, search_start, search_end

    curr_locus = -1
    end_locus = 0

    for i in range(i, len(iter_partitions)):

        pstart, pend = iter_partitions[i]
        f = covar_dict[pstart, pend]

        covar, loci, covar_starts, covar_ends = update_covar_and_loci(covar, loci, end_locus, f)

        autocovars = covar[covar.i == covar.j]

        autocovar = create_fastdict(autocovars.i.values, autocovars.val.values)

        start_locus_index, start_locus = find_start_locus(curr_locus, loci, snp_bottom)
        curr_locus_index, curr_locus = start_locus_index, start_locus
        end_locus_index, end_locus = find_end_locus(i, iter_partitions, loci, snp_last)

        covar_view = covar.values

        v, h = defaultdict(int), defaultdict(int)

        while curr_locus_index < len(loci) and loci.values[curr_locus_index] <= end_locus:

            # precomputed_loci.append(curr_locus)
            curr_locus = loci.values[curr_locus_index]

            if (curr_locus > snp_first or breakpoint_index == 0) and (curr_locus <= snp_last):

                covar_start, covar_end = covar_starts[curr_locus_index], covar_ends[curr_locus_index]

                for x in range(covar_start, covar_end):

                    if covar_view[x][1] <= snp_top:

                        i = covar_view[x][0]
                        j = covar_view[x][1]
                        val = covar_view[x][2]
                        # print("adding", i, j, val)
                        corr_coeff = (val / math.sqrt((autocovar[i] * autocovar[j]))) ** 2

                        v[i] += corr_coeff
                        h[j] += corr_coeff

            curr_locus_index += 1

        precomputed_loci = loci.values[start_locus_index:end_locus_index]

        return v, h, precomputed_loci


if __name__ == "__main__":

    from collections import defaultdict

    assert len(covariance_files) == len(partitions)

    covar_dict = {}
    for f, (start, end) in zip(covariance_files, partitions.itertuples(index=False)):
        covar_dict[start, end] = f

    breakpoints = find_breakpoint_loci(df)
    starts_ends = search_starts_ends(breakpoints, partitions)

    print(starts_ends)

    for breakpoint_index, (search_start, search_end, breakpoint, last_breakpoint, next_breakpoint) in enumerate(starts_ends.itertuples(index=False)):

        result = precompute_values(partitions, search_start, search_end, last_breakpoint, next_breakpoint)
        v, h, precomputed_loci = result
        print("result", precomputed_loci)
        print(len(v))
        print(len(h))


            # print("snp_bottom", snp_bottom)
            # print("snp_top", snp_top)
            # snp_bottom_idx = np.searchsorted(precomputed_loci, snp_bottom, side="right") - 1
            # # print("snp_bottom_idx", snp_bottom_idx)
            # snp_top_idx = np.searchsorted(precomputed_loci, snp_top) - 1
            # print("snp_top_idx", snp_top_idx)


            # zero_metric, loci_to_compute_later, later_bps = compute_zero_metric(precomputed_loci, partitions, breakpoints)
            # metric_sum, nonzero = compute_sum_and_nonzero(loci_to_compute_later, later_bps, covars.i.values, covars.j.values, covars.val.values, autocovar)

            # local_search_right(precomputed_loci, v, h, breakpoint, snp_bottom_idx, snp_top_idx, metric_sum, zero_metric)
            # print("end_locus", end_locus)


        # return v, h
        # print(len(v), len(h))
        # print("here", len(set(v).union(set(h))))
            # if curr_locus


        # for i, (pstart, pend) in enumerate(iter_partitions.itertuples()):

            # to_preread = []

            #     last_pnum = i
            # else:
            #     break
            # for p_num_init in range(0)
            # v, h = precompute_values(covar_dict, snp_top, snp_bottom, snp_first, snp_last)

        # snp_first, snp_last = search_start, search_end


        # snp_bottom, snp_top = break
        # print(last_breakpoint, next_breakpoint)

        # flat.get_final_partitions(self.input_config, self.name, start_search, stop_search)
        # partition i should be straddle snp_first, the latter straddle snp_last
        # temp_partitions = 


        
        # snp_first, snp_last = start_search, stop_search

        # snp_top, snp_bottom = top_bottom_snps()



    # loci = df.pos
    # zero_metric, loci_to_compute_later, later_bps = compute_zero_metric(loci.tolist(), partitions, breakpoints)

    # covars = pd.concat([pd.read_csv(f, sep=" ", usecols=[2, 3, 7], names="i j val".split())
    #                                 for f in covariance_files])

    # autocovar = covars[covars.i.values == covars.j.values].drop("j", 1).set_index("i").squeeze().to_dict()

    # metric_sum, nonzero = compute_sum_and_nonzero(loci_to_compute_later, later_bps, covars.i.values, covars.j.values, covars.val.values, autocovar)

    # assert len(starts_ends) == len(breakpoints), " ".join([str(l) for l in [len(starts_ends),len(breakpoints)]])
    # print(starts_ends)



                        #

                        

            # here we can have teh codez



        # read partition
    # for count, (start, end, breakpoint, next_breakpoint) in enumerate(starts_ends.itertuples(index=False)):

    #     good_covars, good_loci = locs_to_use(covars, count, end, start, next_breakpoint)

    #     if count < len(starts_ends) - 1:
    #         end_locus = partitions.iloc[count + 1, 0]
    #     else:
    #         end_locus = good_loci[np.searchsorted(good_loci, end) - 1]

    #     v, h = initialize_search(good_covars, autocovar)

    #     snp_bottom_idx = np.searchsorted(loci, start, side="right") - 1
    #     # print("snp_bottom_idx", snp_bottom_idx)
    #     snp_top_idx = np.searchsorted(loci, end) - 1
    #     # print("snp_top_idx", snp_top_idx)

    #     local_search_right(good_loci, v, h, breakpoint, snp_bottom_idx, snp_top_idx, metric_sum, zero_metric)
    #     print("end_locus", end_locus)

    # v, h = initialize_search(covars, autocovar)
    # print(len(h))
    
    # breakpoint_loci_local_search = run_local_search_complete(chr_name, breakpoint_loci, begin, end, config, metric_out)

    # print("metric_sum", metric_sum)
    # print("nonzero", nonzero)


    # iterate over covar file to find all needed covariances

    # find autocorrelation for all loci: not that big a dict?
    # then can iterate over each locus to compute and computation is straightforward


# bps_per_loci = next_breakpoint_per_loci(df.pos, breakpoints)

# print(bps_per_loci)
# create array of loci and their breakpoint for simplicity in next algo

# find
# insert = np.searchsorted(df.pos.values, breakpoints, side="right") - 1

