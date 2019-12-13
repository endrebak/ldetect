#!/usr/bin/env python3


import pandas as pd

import numpy as np

# python3 P00_00_partition_chromosome.py <input_genetic_map> <n_individuals_in_ref_panel> <output_file>

# rs187078949 10133 0.0
# rs191522553 10140 0.0
# rs149483862 10286 0.0
# rs150919307 10297 0.0
# rs186644623 10315 0.0
# rs193294418 10345 0.0
# rs185496709 10386 0.0
# rs188771313 10419 0.0
# rs192945962 10425 0.0
# rs184397180 10431 0.0

import sys


# chunk = float(nsnp)/float(N)
# print(chunk)
# chunk = int(math.floor(chunk))
# print(chunk)

# for i in range(chunk):
#     start = i*N
#     end = i*N + N
#     if i == chunk-1:
#         end = len(poss)-1
#         startpos = poss[start]
#         endpos = poss[end]
#         # print >> outfile, startpos, endpos
#         print(str(startpos)+' '+str(endpos), file=outfile)
#         continue
#     startpos = poss[start]
#     endpos = poss[end-1]
#     endgpos = pos2gpos[endpos]
#     test =end+1
#     testpos = poss[test]
#     stop = False
#     while stop == False:
#         if test == len(poss):
#             stop = True
#             continue
#         testpos = poss[test]
#         testgpos = pos2gpos[testpos]
#         df = testgpos-endgpos
#         tmp = math.exp(    -4.0*11418.0*df / (2.0*nind))
#         if tmp < 1.5e-8:
#             stop = True
#         else:
#             test = test+1
#     # print >> outfile, startpos, testpos
#     print(str(startpos)+' '+str(testpos), file=outfile)




def partition_chromosomes(genetic_map, n_individuals):

    df = pd.read_csv(genetic_map, sep=" ", header=None, names="pos gpos".split(), nrows=None, usecols=[1, 2])
    chunksize = 5000

    print(df)

    df = df.drop_duplicates("pos")

    print(df)

    n_chunks = int(np.floor((len(df)/float(chunksize))))
    print(n_chunks)

    # print(chunk)

    print(np.array_split(df, n_chunks))

if __name__ == "__main__":
    partition_chromosomes(sys.argv[1], sys.argv[2])

# infile = gzip.open(sys.argv[1])
# nind = float(sys.argv[2])
# outfile = open(sys.argv[3], "w")
# N = 5000


# pos2gpos = dict()
# poss = list()
# line = infile.readline()
# while line:
#     line = line.strip().split()
#     pos = int(line[1])
#     gpos = float(line[2])
#     pos2gpos[pos] = gpos
#     poss.append(pos)
#     line = infile.readline()

# print(len(poss))
# nsnp = len(poss)

