#!/usr/bin/env python3.4

import sys
import matplotlib.pylab as pylab
import matplotlib.pyplot as pt
import numpy as np
import math 
import ldetect.baselib.flat_file_consts as cnst
import ldetect.baselib.flat_file as flat
import ldetect.baselib.filters as filt
import ldetect.baselib.read_data as rd

# GRAPH ALL:

def plot_num_of_minima_vs_filter_width(graphs):
	data = [(el['width'],len(el['filtered_minima_ind'])) for el in graphs]
	data = sorted(data, key=lambda x: x[0])

	x = [2*el[0]+1 for el in data]
	y = [el[1] for el in data]

	pt.figure()

	pt.xlabel('filter width')
	pt.ylabel('# of minima')
	pt.title('# of minima vs filter width')

	pt.plot(x,y)
	# pt.show()

def plot_specgram(in_array):
	pt.figure()
	pylab.specgram(in_array, NFFT=256, noverlap=224)
	
	pt.xlabel('SNP #')
	pt.ylabel('f')
	pt.suptitle('Spectrogram')
	pt.title('Color is |A|^2')

def plot_spectrum(in_array, plot):
	w_param = 4
	
	freqs = np.fft.fftfreq(len(in_array)*w_param, 1)
	idx = np.argsort(freqs)
	ps = np.abs(np.fft.fft(in_array,len(in_array)*w_param))**2
	
	plot.set_yscale('log')
	plot.set_xlabel('f')
	plot.set_ylabel('|A|^2')
	plot.set_title('Frequency Spectrum')

	plot.plot(freqs[idx], ps[idx])

def plot_spectrum_single_fig(in_array, plot):
	w_param = 4
	
	freqs = np.fft.fftfreq(len(in_array)*w_param, 1)
	idx = np.argsort(freqs)
	ps = np.abs(np.fft.fft(in_array,len(in_array)*w_param))**2
	plot.plot(freqs[idx], ps[idx])

	plot.yscale('log')
	plot.xlabel('f')
	plot.ylabel('|A|^2')
	plot.title('Frequency Spectrum')

def add_plot(g, np_init_array, x, axarr, i, show_minima=True):
	# pt.plot(x,g['filtered'])
	limit_val = np.amax(g['filtered'])

	axarr[i,0].set_xlabel('SNP #')
	axarr[i,0].set_ylabel('Sum')
	axarr[i,0].set_title('Sum Vector')

	axarr[i,0].plot(x,g['filtered'])

	if show_minima:
		# Draw the minima as borders
		minima_loc = [0]*len(np_init_array)
		for ind in g['filtered_minima_ind']:
			minima_loc[ind] = limit_val

		axarr[i,0].plot(x,minima_loc)	

		# Draw borders when it's a uniform distribution
		uniform_loc = [0]*len(np_init_array)
		cnt = 0
		for ind in range(0,len(np_init_array), int(len(np_init_array)/(len(g['filtered_minima_ind'])+1))):
			uniform_loc[ind] = limit_val/10
			cnt+=1

		axarr[i,0].plot(x,uniform_loc)	


	axarr[i,1].set_xlabel('Distance between minima')
	axarr[i,1].set_ylabel('# of minima')
	axarr[i,1].set_title('Histogram')
	axarr[i,1].hist(g['filtered_minima_vals'], 40)

	plot_spectrum(g['window'], axarr[i,2])
	plot_spectrum(g['window_moving_avg'], axarr[i,2])


def plot_filtered_and_minima_distribution_and_filter_spectrum(np_init_array, graphs, np_init_array_x=None, draw_raw_data=True):
	if np_init_array_x is None:
		x = list(range(0,len(np_init_array))) # done here just so we don't do it every time a graph is added
	else:
		x = np_init_array_x 

	extra_graph = 0
	if draw_raw_data:
		extra_graph = 1

	f, axarr = pt.subplots(len(graphs)+extra_graph, 3, sharex='col') #, sharey='col')

	if draw_raw_data:
		graphs_dummy = filt.apply_filters(np_init_array, 0, 0, 1)
		add_plot(graphs_dummy[0], np_init_array, x, axarr, 0, False)

	for i,g in enumerate(graphs):
		add_plot(g, np_init_array, x, axarr, i+extra_graph)

def plot_all(np_init_array, graphs, np_init_array_x=None):
	print(len(graphs))

	plot_specgram(np_init_array)
	plot_num_of_minima_vs_filter_width(graphs)
	plot_filtered_and_minima_distribution_and_filter_spectrum(np_init_array, graphs, np_init_array_x)

	pt.show()

def standard_run(np_init_array, np_init_array_x, start, stop, step):
	# Interactive plots
	graphs = filt.apply_filters(np_init_array, start, stop, step)
	# graphs = filt.apply_filters(np_init_array, int(sys.argv[1]), int(sys.argv[2])+1, int(sys.argv[3]))
	
	for g in graphs:
		flat.print_log_msg('indices'+repr(g['width'])+repr(g['filtered_minima_ind']))
		loci = filt.get_minima_loc(g, np_init_array_x)
		flat.print_log_msg('loci'+repr(g['width'])+repr(loci))
		
	plot_all(np_init_array, graphs, np_init_array_x)

	# pt.figure()
	# plot_spectrum_single_fig(np_init_array, pt)
	# pt.figure()
	# plot_spectrum_single_fig(graphs[0]['filtered'], pt)
	# pt.show() 

def main():
	'''
	Run with at least:
	a) 1 argument, defining the central filter width of analysis area
	b) 3 arguments, defining the start, stop, and step of series of filter widths

	optional last argument is file name of data to read into memory (if it hasn't alreay been read)
	'''
	if sys.argv[1] == 'help':
		print(main.__doc__)
		return
	# max_w = 10000
	# vals = []
	# for width in range(1, max_w):
	# 	vals.append(filt.apply_filter_get_minima(np_init_array, width))
	# 	print(width)

	# np_temp_array = np.array(vals)

	# minima = sig.argrelextrema(np_temp_array, np.greater)[0]

	# print(minima)

	global np_init_array, np_init_array_x

	try:
		np_init_array
		np_init_array_x
	except NameError:
		flat.print_log_msg('Reading data')
		np_init_array, np_init_array_x = rd.read_data(sys.argv[len(sys.argv)-1])

	relative_width = 0.5

	if len(sys.argv) > 3:
		start, stop, step = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
	else:
		center_val = int(sys.argv[1])
		start, stop, step = math.floor(center_val-relative_width*center_val), math.ceil(center_val+relative_width*center_val), math.floor(2*relative_width*center_val/6)
	
	standard_run(np_init_array, np_init_array_x, start, stop, step)

	flat.print_log_msg('Done')

if __name__ == '__main__':
	main()


# freqs = np.fft.fftfreq(len(g['filtered'])*2, 1)
# idx = np.argsort(freqs)
# ps = np.abs(np.fft.fft(a,len(g['filtered'])*2))**2
# axarr[i,3].plot(freqs[idx], ps[idx])
# axarr[i,3].set_yscale('log')

# axarr[i,0].plot(x,g['filtered'])
# axarr[i,1].hist(g['filtered_minima_ind_diff'],bins=1000)

# --------------------------------------------------------

# axarr[len(graphs)].plot(x,np_init_array)

# axarr[len(graphs),0].plot(x,np_init_array)
# axarr[len(graphs),1].hist(g['filtered_minima_ind_diff'],bins=1000)

# minima_a = sig.argrelextrema(np_init_array, np.less)[0]
# minima_a_diff = [minima_a[i]-minima_a[i-1] for i in range(1, len(minima_a))]
# axarr[len(graphs),1].hist(minima_a_diff,bins=1000)

# pt.plot(x,np_init_array,'.',label='orig')

# --------------

# pt.show()

# --------------

# testing Smooth w/ gaussian:
# width = 100
# b=sig.gaussian(2*width+1, width/3); ga=filters.convolve1d(np_init_array, b/b.sum()); pt.plot(ga); pt.show()

# testing Smooth w/ hamming:
# In [232]: a = sig.firwin(200, cutoff=0.001)
# In [233]: ga=filters.convolve1d(init_array_np, a)

# testing SPECTROGRAM - tweak NFFT and noverlap
# specgram = mlab.specgram(in_list, NFFT=128, noverlap=112); pt.pcolormesh(specgram[0]); pt.colorbar(); pt.show()

# testing CWT:
# np_init_array = np.array(in_list)
# wavelet = sig.ricker
# widths = np.arange(1, 1000)
# cwtmatr = sig.cwt(np_init_array, wavelet, widths); pt.pcolormesh(cwtmatr); pt.colorbar(); pt.show()

# ---------------