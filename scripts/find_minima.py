import pandas as pd

from sys import argv
import scipy.ndimage.filters as filters
import scipy.signal as sig

import math
sqrt = math.sqrt
import numpy as np

from time import time

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

    # print("loci", loci)
    # print("bps", bps)
    # print("covars", covars)

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












if __name__ == "__main__":

    breakpoints = find_breakpoint_loci(df)

    zero_metric, loci_to_compute_later, later_bps = compute_zero_metric(df.pos.tolist(), partitions, breakpoints)

    covars = pd.concat([pd.read_csv(f, sep=" ", usecols=[2, 3, 7], names="i j val".split())
                                    for f in covariance_files])

    autocovar = covars[covars.i.values == covars.j.values].drop("j", 1).set_index("i").squeeze().to_dict()

    # print(autocovar)

    metric_sum, nonzero = compute_sum_and_nonzero(loci_to_compute_later, later_bps, covars.i.values, covars.j.values, covars.val.values, autocovar)

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

