import ldetect.baselib.flat_file_consts as cnst
import ldetect.baselib.flat_file as flat

import gzip
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as pt
import csv

def main():
	name = 'chr1'

	flat.print_log_msg('Starting run')
	x, y, pairs = flat.read_hotspots(cnst.const['genetic_maps']['root']+cnst.const['genetic_maps']['file_base']+name+cnst.const['genetic_maps']['ext'])
	
	flat.print_log_msg('Plotting')
	pt.plot(x,y)
	fig = pt.gcf()
	fig.set_size_inches((40,30))

	pt.xlabel('SNP #')
	pt.ylabel('Hotspot val')
	pt.title('Hotspots')

	pt.savefig('genetic_maps_output.png')

	pt.clf()
	pt.plot(x[5000:10000],y[5000:10000])
	fig = pt.gcf()
	fig.set_size_inches((40,30))

	pt.xlabel('SNP #')
	pt.ylabel('Hotspot val')
	pt.title('Hotspots zoomed')

	pt.savefig('genetic_maps_zoomed_output.png')

	flat.print_log_msg('Done')
	# pt.show

if __name__ == '__main__':
	main()