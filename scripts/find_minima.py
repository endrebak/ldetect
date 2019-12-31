import pandas as pd

from sys import argv
import scipy.ndimage.filters as filters
import scipy.signal as sig

import math
import numpy as np

from time import time

f = argv[1]

df = pd.read_table(f, header=None, names="pos val".split())

# def _apply_filter(a, width):

#     val = (2 * width + 1)
#     moving_avg_a = np.ones(val) * 1/(val)
#     # start = time()
#     a = sig.get_window('hanning', val)
#     # end = time()
#     # print("took", end - start)
#     return a

def _apply_filter(arr, width):

    val = (2 * width + 1)
    moving_avg_a = np.ones(val) * 1/(val)
    # start = time()
    a = sig.get_window('hanning', val)

    ga = filters.convolve1d(arr, a/a.sum())

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


def _width_raw(a, n_bpoints):

    "When looking up a index, you use that index as the width for the filter"

    "binary search but apply the filter for each time"

    search_val = n_bpoints
    # print("search_val", search_val)

    filtered = _apply_filter(a, search_val)
    idx = np.searchsorted(filtered, search_val, side="left") - 1

    return len(filtered) - idx


def _trackback(array, search_val, start_search):

    step_coarse = 200
    delta_coarse = 20
    for i in range(start_search + step_coarse, start_search + delta_coarse, step_coarse):
        pass
        



def find_minima(df, n_snps=50):

    print(len(df))
    print(n_snps - 1)

    n_bpoints = int(math.ceil(len(df) / n_snps) - 1)

    print("Number of breakpoints", n_bpoints)

    search_value = 1000
    pos = df.pos
    val = df.val
    end = _find_end(val, search_value, n_bpoints)

    print("End", end)

    found_width_raw = _width_raw(val, n_bpoints)
    print("found_width_raw", found_width_raw)
    found_width = end - found_width_raw
    print("found_width:", found_width)

    found_width_trackback_raw = _trackback(val, search_value, found_width_raw)



find_minima(df)
