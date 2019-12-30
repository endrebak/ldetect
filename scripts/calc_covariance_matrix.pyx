#!/usr/bin/env python3
#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, 

import sys, math, gzip
import numpy as np
# cimport numpy as cnp
import pandas as pd

from time import time

from libc.math cimport exp, fabs
from libc.stdint cimport int8_t

# # calculate Wen/Stephens shrinkage LD estimate
# gmapfile = sys.argv[1] # genetic map
# indfile = sys.argv[2] #list of individuals
# # NE = 11418.0
# # CUTOFF = 1e-7
# outfile = sys.argv[5] # outfile file

# NE = float(sys.argv[3])
# CUTOFF = float(sys.argv[4])


cpdef calc_covar(gmapfile, indfile, double NE, double CUTOFF):

    inds = pd.read_table(indfile, header=None, squeeze=True).to_list()

    haps= list()
    nind = len(inds)
    s = 0

    for _i in range(1, 2*nind):
            s = s+ 1.0/float(_i)
    nind = float(nind)

    s = 1/s
    #print "s", s


    cdef:
        int i, j, k, pos1, pos2, len_allpos, toofar, n11, n10, n01, len_g1_int
        double gpos1, gpos2, ee_const, len_g1, ee, df, theta, f11, f1, f2, b, Ds, Ds2, D

    pos2gpos = pd.read_table(gmapfile, index_col=0, header=None, sep=" ", squeeze=True)
    pos2gpos = pos2gpos

    df = pd.read_table(sys.stdin, header=None)
    allpos = df.pop(0) #.tolist()
    allrs = df.pop(1) #.tolist()
    # haps = df.astype(np.int8)
    haps = df.values

    pos2gpos = pos2gpos[allpos]

    records = []

    len_g1_int = haps.shape[1]
    len_g1 = float(haps.shape[1])
    ee_const = NE*4.0 / (2.0*nind)
    len_allpos = len(allpos)
    theta = s/(2.0*float(nind)+s)

    # g1_range = range(len_g1)

    cdef double[:] pos2gpos_view
    cdef int[:] allpos_view
    cdef int[:] g1_view, g2_view
    cdef int[:, :] haps_view
    # cdef cnp.ndarray[int8_t, ndim=2] haps_arr
    # haps_arr = haps
    # cdef int[::] haps_view
    # g1_range_view
    pos2gpos_view = pos2gpos
    allpos_view = allpos
    haps_view = haps.values

    for i in range(len_allpos):

            # print("  i", i)

            pos1 = allpos_view[i]

            gpos1 = pos2gpos_view[i]

            toofar = 0
            j = i

            while j < len_allpos and not toofar:
                    pos2 = allpos_view[j]
                    gpos2 = pos2gpos_view[j]

                    df = gpos2 - gpos1
                    ee = exp( - df * ee_const)

                    if ee < CUTOFF:
                            toofar = 1
                            j = j+1
                            continue

                    # g1_view = haps_view[i]

                    # g2_view = haps_view[j]

                    # prin=(g1_view)

                    n11, n01, n10 = 0, 0, 0
                    for k in range(len_g1_int):
                        if haps_view[i][k] == 1 and haps_view[j][k] == 1:
                            n11 += 1
                        elif haps_view[i][k] == 0 and haps_view[j][k] == 1:
                            n01 += 1
                        elif haps_view[i][k] == 1 and haps_view[j][k] == 0:
                            n10 += 1

                    # haps_compare = pd.concat([g1, g2], axis=1)
                    # haps_compare.columns = [0, 1]
                    # haps_compare = haps_compare.groupby([0, 1]).size().astype(float).to_dict()
                    # n11 = haps_compare.get((1, 1), 0)
                    # n10 = haps_compare.get((1, 0), 0)
                    # n01 = haps_compare.get((0, 1), 0)

                    f11 = n11/len_g1
                    f1 = (n11+n10)/len_g1
                    f2 = (n11+n01)/len_g1
                    D = f11 - f1*f2
                    Ds = D*ee
                    Ds2 = (1-theta)*(1-theta)*Ds

                    if fabs(Ds2) < CUTOFF:
                            j = j+1
                            continue
                    if i == j:
                            Ds2 = Ds2 + (theta/2.0)*(1-theta/2.0)

                    result = (allrs[i], allrs[j], pos1, pos2, gpos1, gpos2, D, Ds2)
                    records.append(result)
                    j = j+1

    df = pd.DataFrame.from_records(records)
    return df

