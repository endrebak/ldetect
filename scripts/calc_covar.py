#!/usr/bin/env python3

import sys
import pandas as pd
# import numpy as np
# import pandas as pd

# calculate Wen/Stephens shrinkage LD estimate
gmapfile = sys.argv[1] # genetic map
indfile = sys.argv[2] #list of individuals
# NE = 11418.0
# CUTOFF = 1e-7
# outfile = sys.argv[5] # outfile file

NE = float(sys.argv[3])
CUTOFF = float(sys.argv[4])


from ldetect2.src.calc_covar import calc_covar



df = pd.read_table(sys.stdin, header=None)
calc_covar(df, gmapfile, indfile, NE, CUTOFF)


# print("writing to", outfile)
# df.to_csv(outfile, sep="\t", index=False, header=False)

