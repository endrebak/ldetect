# Pickle input as alternative

import matplotlib.pyplot as pt
import pickle

def generate_data(max_w):
	# Actually calculate the values (or input from file)
	vals = []
	for width in range(1, max_w):
		val = filt.apply_filter_get_minima(np_init_array, width)
		vals.append(val)
		print(width, val)

	return vals

def export_data(vals, out_fname):
	with open(out_fname, 'wb') as f_out:
		pickle.dump(vals, f_out)

def import_data(in_fname):
	with open(in_fname, 'rb') as f_in:
		vals = pickle.load(f_in)

	return vals

def analyze_data(vals):
	# Find where jumps are
	jump_locs = [i for i in range(1,len(vals)) if vals[i] > vals[i-1]]

	# Find where they start (trace back)
	jump_pairs = []

	for j in jump_locs:
		val = vals[j]
		for p in range(j-1, 0, -1):
			if vals[p] >= val:
				jump_pairs.append(p)
				break

	# How long are they?
	jump_lengths = []

	for i in range(0, len(jump_pairs)):
		jump_lengths.append(jump_locs[i] - jump_pairs[i])

	return jump_locs, jump_pairs, jump_lengths

def plot_analysis_and_data(jump_locs, jump_pairs, jump_lengths, vals):
	# Output & plot
	print('jump_lengths: ', jump_lengths)
	print('max(jump_lengths):', max(jump_lengths))

	# For legibility plot all at same height (500)
	y1 = [500]*len(jump_locs)
	y2 = [500]*len(jump_pairs)

	# Also plot near to the curve, so it can be seen when zooming in
	y3 = [1.02*vals[i] for i in jump_locs]
	y4 = [1.02*vals[i] for i in jump_pairs]

	pt.plot(vals)
	pt.scatter(jump_locs, y1, c='r')
	pt.scatter(jump_pairs, y2, c='b')
	pt.scatter(jump_locs, y3, c='r')
	pt.scatter(jump_pairs, y4, c='b')

	pt.xlabel('filter width')
	pt.ylabel('# of minima')
	pt.title('# of minima vs filter width')

	pt.show()

if __name__ == '__main__':
	# do the steps:
	action = 'import'
	in_fname = 'num_of_minima_1.pickle'
	out_fname = 'num_of_minima_new.pickle'
	max_w = 10000

	if action == 'import':
		# import
		vals = import_data(in_fname)
	elif action == 'generate':
		# generate
		vals = generate_data(max_w)
		# export
		export_data(vals, out_fname)
	else:
		raise Exception('Error: unknown action!')

	# Analyze
	jump_locs, jump_pairs, jump_lengths = analyze_data(vals)

	# output+plot
	plot_analysis_and_data(jump_locs, jump_pairs, jump_lengths, vals)

