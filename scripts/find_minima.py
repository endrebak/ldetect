import pandas as pd

from sys import argv
import scipy.ndimage.filters as filters
import scipy.signal as sig

import math
sqrt = math.sqrt
import numpy as np

from time import time, sleep


f = argv[1]

df = pd.read_table(f, header=None, names="pos val".split())

f2 = argv[2]

partitions = pd.read_table(f2, header=None, sep=" ")

covariance_files = argv[3:]

# def _apply_filter(a, width):

#     val = (2 * width + 1)
#     moving_avg_a = np.ones(val) * 1/(val)
#     # start = time()
#     a = sig.get_window('hanning', val)
#     # end = time()
#     # print("took", end - start)
#     return a

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
    df = pd.concat([search_starts, search_ends, bps, bps_next], axis=1)

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

    
if __name__ == "__main__":

    breakpoints = find_breakpoint_loci(df)

    loci = df.pos
    zero_metric, loci_to_compute_later, later_bps = compute_zero_metric(loci.tolist(), partitions, breakpoints)

    covars = pd.concat([pd.read_csv(f, sep=" ", usecols=[2, 3, 7], names="i j val".split())
                                    for f in covariance_files])

    autocovar = covars[covars.i.values == covars.j.values].drop("j", 1).set_index("i").squeeze().to_dict()

    metric_sum, nonzero = compute_sum_and_nonzero(loci_to_compute_later, later_bps, covars.i.values, covars.j.values, covars.val.values, autocovar)

    starts_ends = search_starts_ends(breakpoints, partitions)
    assert len(starts_ends) == len(breakpoints), " ".join([str(l) for l in [len(starts_ends),len(breakpoints)]])
    print(starts_ends)

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

    def read_partition(p_num):

        pass

    def find_start_locus(curr_locus, loci, snp_bottom):

        if curr_locus < 0:
            for i, locus in enumerate(loci):
                if locus >= snp_bottom:
                    return i, locus
        else:
            return loci.index(curr_locus), curr_locus

    def find_end_locus(partitions, loci, snp_last):

        if p_num + 1 < len(partitions):
            end_locus = partitions.iloc[p_num + 1, 0] # start of next partition
            end_locus_index = -1
            return end_locus_index, end_locus
        else:
            for i in reversed(range(0, len(loci))):
                if loci[i] <= snp_last:
                    return i, end_locus

    def update_covar_and_loci(covar, loci, pos_to_delete_below, partition_to_add):
        pass

def read_partitions(partitions):
    snp_bottom = 0
    snp_last = 0
    for p_num, partition in enumerate(partitions[:-1]):

        to_preread = []
        last_pnum = -1
        if snp_bottom >= partition[0]:
            to_preread.append(p_num)
            last_pnum = p_num
        else:
            break

        covar, loci = read_partition(to_preread)

        curr_locus = -1
        for p_num, partition in enumerate(partitions[last_pnum + 1:], last_pnum + 1):

            covar, loci = update_covar_and_loci(covar, loci, pos_to_delete_below, partition_to_add)

            curr_locus_index, curr_locus = find_start_locus(curr_locus, loci, snp_bottom)
            end_locus_index, end_locus = find_end_locus(curr_locus, loci, snp_last)

            # here we can have teh codez
            while curr_locus <= end_locus:
                loci_



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

