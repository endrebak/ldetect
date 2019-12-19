#!/usr/bin/env python3

import sys, math
import pandas as pd

infile = snakemake.input["genetic_map"]
nind = int(open(snakemake.input["number_individuals"]).readline().strip())
outfile = snakemake.output[0]

N = 5000

indf = pd.read_table(infile, sep=" ", header=None)
# print(indf)
poss = indf[0].tolist()
pos2gpos = {p: g for p, g in zip(poss, indf[1].tolist())}

nsnp = len(poss)

chunk = float(nsnp)/float(N)
chunk = int(math.floor(chunk))

results = []

for i in range(chunk):
    start = i*N
    end = i*N + N
    if i == chunk-1:
        end = len(poss)-1
        startpos = poss[start]
        endpos = poss[end]
        results.append((startpos, endpos))
        continue
    startpos = poss[start]
    endpos = poss[end-1]
    endgpos = pos2gpos[endpos]
    test =end+1
    testpos = poss[test]
    stop = False
    while stop == False:
        if test == len(poss):
            stop = True
            continue
        testpos = poss[test]
        testgpos = pos2gpos[testpos]
        df = testgpos-endgpos
        tmp = math.exp(    -4.0*11418.0*df / (2.0*nind))
        if tmp < 1.5e-8:
            stop = True
        else:
            test = test+1

    results.append((startpos, testpos))

df = pd.DataFrame.from_records(results)
df.to_csv(outfile, sep=" ", index=False, header=False)
