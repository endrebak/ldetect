import numpy as np
import csv
import sys
import gzip
import os.path

def read_data(fname):
	'''
	Uses read_data_raw and converts to numpy arrays
	'''
	in_list, in_list_x = read_data_raw(fname)

	np_init_array_x = np.array(in_list_x)
	np_init_array = np.array(in_list)

	return np_init_array, np_init_array_x

def read_data_np(fname):
	'''
	This is here just for naming consistency - equivalent to calling read_data(fname)
	'''
	return read_data(fname)

def read_data_raw(fname):
	'''
	Reads data from fname (gzipped or plaintext) and returns raw python list 
	'''

	extension = os.path.splitext(fname)[1]

	if extension in ('.gz','.GZ','.Gz','.gZ'):
		return read_data_raw_zipped(fname) 
	else:
		return read_data_raw_plaintext(fname)

def read_data_raw_plaintext(fname):	
	with open(fname, 'r') as f_in:
		return csv_read(f_in)

	return None, None # If error w/ file

def read_data_raw_zipped(fname):
	with gzip.open(fname, 'rt') as f_in:
		return csv_read(f_in)

	return None, None # If error w/ file

def csv_read(f_in):
	csv_reader=csv.reader(f_in, delimiter='\t')
	in_list = []
	in_list_x = []
	for row in csv_reader:
		in_list_x.append(int(row[0]))
		in_list.append(float(row[1]))

	return in_list, in_list_x	

if __name__ == '__main__':
	print(sys.argv[1])
	np_init_array, np_init_array_x = read_data(sys.argv[1])