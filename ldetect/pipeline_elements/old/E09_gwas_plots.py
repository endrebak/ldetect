#!/usr/bin/env python3.4

import matplotlib.pyplot as pt
import csv
import numpy as np
import pickle
import itertools

import commanderline.commander_line as cl
import ldetect.baselib.flat_file_consts as cnst
import ldetect.baselib.flat_file as flat

def plot_gwas(input_bed_fname, input_pickle_fname):
	with open(input_bed_fname, 'rt') as f_in:
		csv_reader = csv.reader(f_in, delimiter='\t')
		x = []
		y = []
		for row in csv_reader:
			x.append(int(row[1]))
			y.append(float(row[4]))

		np_array_x = np.array(x)
		np_array_y = np.array(y)

		colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y'])
		markers = itertools.cycle(['o', 'x', '+', '*', '<', '>'])

		pt.scatter(np_array_x, np_array_y, c=next(colors), label='p-values', marker=next(markers))

		with open(input_pickle_fname, 'rb') as f_in:
			loaded = pickle.load(f_in)

			# print(loaded)

			dist = -0.1

			for key in loaded:
				# print()
				# print(key, loaded[key])
				try:
					if 'loci' in loaded[key]:
						dist -= 0.05
						breakpoint_list_y = [dist for a in loaded[key]['loci']]
						pt.scatter(np.array(loaded[key]['loci']), np.array(breakpoint_list_y), c=next(colors), label=key, marker=next(markers))
				except TypeError as e:
					print(e)

		pt.legend()
		pt.show()

def plot_gwas_w_default_paths(chr_name, input_pickle_fname):
	plot_gwas(cnst.const['gwas_plots']['height']['root']+cnst.const['gwas_plots']['height']['file_prefix']+chr_name+cnst.const['gwas_plots']['height']['file_suffix'], input_pickle_fname)


if __name__ == '__main__':
    cl.commander_line((plot_gwas_w_default_paths, plot_gwas))