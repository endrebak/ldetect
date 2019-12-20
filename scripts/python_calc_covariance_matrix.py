#!/usr/bin/env python3

import sys, math, gzip
import numpy as np
import pandas as pd

from time import time

# calculate Wen/Stephens shrinkage LD estimate
gmapfile = sys.argv[1] # genetic map
indfile = sys.argv[2] #list of individuals
# NE = 11418.0
NE = float(sys.argv[3])
# CUTOFF = 1e-7
CUTOFF = float(sys.argv[4])
outfile = sys.argv[5] # outfile file


inds = pd.read_table(indfile, header=None, squeeze=True).to_list()

haps= list()
theta = 0
nind = len(inds)
s = 0
for i in range(1, 2*nind):
        s = s+ 1.0/float(i)
nind = float(nind)

s = 1/s
#print "s", s
theta = s/(2.0*float(nind)+s)
print(theta)


pos2gpos = pd.read_table(gmapfile, index_col=0, header=None, sep=" ", squeeze=True)
# pos2gpos = pos2gpos.to_dict()



print(list(pos2gpos.items())[:5])
print(len(pos2gpos))

df = pd.read_table(sys.stdin, header=None)
allpos = df.pop(0) # .tolist()
allrs = df.pop(1)  # .tolist()
haps = df.astype(np.int8)

pos2gpos = pos2gpos[allpos]


records = []
len_g1 = float(haps.shape[1])
ee_const = NE * 4.0 / (2.0 * nind)
for i in range(len(allpos)):

        # if i == 1:
        #         raise Exception("Hi")
        # print("-----")

        # print("i", i, i/len(allpos), len(allpos))
        # print(time())

        pos1 = allpos[i]

        # print(pos1)

        gpos1 = pos2gpos[pos1]

        # print("  gpos1", gpos1)

        toofar = False
        j = i
        # print("  len(allpos)", len(allpos))
        # print("j, len(allpos)")
        # print(j, len(allpos))
        while j < len(allpos) and toofar == False:
                # whole_start = time()
                # print("  i", i, "j", j)
                pos2 = allpos[j]
                gpos2 = pos2gpos[pos2]

                # print("  gpos2", gpos2)
                df = gpos2-gpos1
                # print("  df", df)
                # print("  NE", NE)
                # print("  inds", len(inds))
                ee = math.exp( - df * ee_const)
                # print("  ee", ee)
                # print("  CUTOFF", CUTOFF)
                if ee < CUTOFF:
                        toofar = True
                        j = j+1
                        continue
                g1 = haps.iloc[i]
                # print("  g1", g1)
                g2 = haps.iloc[j]
                haps_compare = pd.concat([g1, g2], axis=1)
                haps_compare.columns = [0, 1]
                haps_compare = haps_compare.groupby([0, 1]).size().astype(float).to_dict()
                n11 = haps_compare.get((1, 1), 0)
                n10 = haps_compare.get((1, 0), 0)
                n01 = haps_compare.get((0, 1), 0)
                # end = time()
                # print("took", end - start)
                # print("  g1", g1)
                # for k in range(len(g1)):
                #         # print("  k", k)
                #         if g1[k] == "1" and g2[k] == "1":
                #                 n11 = n11+1
                #         elif g1[k] == "0" and g2[k] == "1":
                #                 n01 = n01 +1
                #         elif g1[k] == "1" and g2[k] == "0":
                #                 n10 = n10 +1
                f11 = n11/len_g1
                f1 = (n11+n10)/len_g1
                f2 = (n11+n01)/len_g1
                D = f11 - f1*f2
                Ds = D*ee
                Ds2 = (1-theta)*(1-theta)*Ds

                # whole_end = time()
                # print("whole took", whole_end - whole_start)
                if math.fabs(Ds2) < CUTOFF:
                        j = j+1
                        continue
                if i == j:
                        Ds2 = Ds2 + (theta/2.0)*(1-theta/2.0)

                result = (allrs[i], allrs[j], pos1, pos2, gpos1, gpos2, D, Ds2)
                print(result)
                records.append(result)
                # print("  j", j)
                j = j+1

df = pd.DataFrame.from_records(records)

df.to_csv(outfile, sep="\t", index=False, header=False)
