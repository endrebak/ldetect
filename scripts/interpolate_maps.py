#!/usr/bin/python

import pandas as pd
import sys

try:
    bedfile = snakemake.input.bed #.bed
    mapfile = snakemake.input.mapfile #input map file, either the HapMap map or the 1000 Genomes OMNI map
    outfile = snakemake.output[0] #output style: [rs] [pos] [genetic pos]
except NameError:
    bedfile = sys.argv[1] #.bed
    mapfile = sys.argv[2] #input map file, either the HapMap map or the 1000 Genomes OMNI map
    outfile = sys.argv[3] #output style: [rs] [pos] [genetic pos]

bdf = pd.read_table(bedfile, header=None)
# print(bdf)
# rsin = bdf[3].tolist()
posin = bdf[1].tolist()

mdf = pd.read_table(mapfile, header=None, skiprows=1)
# print(mdf)

mappos = mdf[1].tolist()
mapgpos = mdf[3].tolist()

results = []

index1 = 0
index2 = 0

while index1 < len(posin):
    pos = posin[index1]
    # rs = rsin[index1]
    if pos == mappos[index2]:
        #the 1000 Genomes site was genotyped as part of the map
        results.append((pos, mapgpos[index2]))
        index1 = index1 + 1
    elif pos < mappos[index2]:
        #current position in interpolation before marker
        if index2 == 0:
            #before the first site in the map (genetic position = 0)
            # print(rs, pos, mapgpos[index2], file=outfile)
            results.append((pos, mapgpos[index2]))
            index1 = index1 + 1
        else:
            #interpolate
            prevg = mapgpos[index2 - 1]
            prevpos = mappos[index2]
            frac = (float(pos) - float(mappos[index2 - 1]))/ (float(mappos[index2]) - float(mappos[index2 - 1]))
            tmpg = prevg + frac * (mapgpos[index2] - prevg)
            results.append((pos, tmpg))
            index1 = index1 + 1
    elif pos > mappos[index2]:
        #current position in interpolation after marker
        if index2 == len(mappos) - 1:
            #after the last site in the map (genetic position = maximum in map,
            # note could try to extrapolate based on rate instead)
            results.append((pos, mapgpos[index2]))
            index1 = index1 + 1
        else:
            index2 = index2 + 1

df = pd.DataFrame.from_records(results)
df.to_csv(outfile, sep=" ", header=False, index=False)

