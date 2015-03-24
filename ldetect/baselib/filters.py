import numpy as np
import scipy.signal as sig
import scipy.ndimage.filters as filters

def apply_filter_get_minima(np_init_array, width):
	el = apply_filter(np_init_array, width)

	return len(el['filtered_minima_ind'])

def apply_filter_get_minima_ind(np_init_array, width):
	el = apply_filter(np_init_array, width)

	return el['filtered_minima_ind']

def get_minima_loc(g, np_init_array_x):
	loci = [np_init_array_x[l] for l in g['filtered_minima_ind']]
	
	return loci

def apply_filter(np_init_array, width):
	# a=sig.gaussian(2*width+1, width/3)
	# a=sig.firwin(width, cutoff=0.001)

	# a=sig.get_window(('gaussian',width/4),2*width+1)
	# a=sig.get_window('hamming',2*width+1)
	moving_avg_a=np.array([1/(2*width+1)]*(2*width+1))
	a=sig.get_window('hanning',2*width+1)
	# a = moving_avg_a
	
	ga=filters.convolve1d(np_init_array, a/a.sum())
	# ga = sig.convolve(np_init_array, a/a.sum(), mode='valid')

	minima_a = sig.argrelextrema(ga, np.less)[0]

	print(2*width+1, len(minima_a))

	minima_a_vals = [ga[i] for i in minima_a]

	el = {  'width':width, 
			'window':a, 
			'window_moving_avg':moving_avg_a,
			'filtered':ga, 
			'filtered_minima_ind':minima_a, 
			'filtered_minima_vals':minima_a_vals, 
	}

	return el

def apply_filters(np_init_array, first, last, step):
	graphs = []
	for width in range(first, last+1, step):
		el = apply_filter(np_init_array, width)

		graphs.append(el)
		#pt.plot(x,ga,'-',label=str(width))
		# print(width)

	return graphs