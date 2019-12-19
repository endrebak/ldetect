# from . import flat_file_consts as cnst

import sys
import csv
import gzip
import time
import math
import bisect

def get_final_partitions(input_config, name, snp_first, snp_last):
	if snp_first>snp_last:
		raise Exception('Error: snp_first is larger than snp_last!')

	print_log_msg('Reading partitions file')
	partitions = read_partitions(name, input_config)

	print_log_msg('Getting relevant partitions')
	partitions = relevant_subpartitions(partitions, snp_first, snp_last)

	print(partitions)

	if len(partitions)<=0:
		raise Exception('Error: There are no relevant subpartitions')

	return partitions

def read_partitions(name, input_config):
	partitions = []
	with open(input_config['partitions_dir']+name+input_config['partitions_file_suffix']) as f_part:
		csv_reader = csv.reader(f_part, delimiter=input_config['partitions_delim'])
		for row in csv_reader:
			partitions.append((int(row[input_config['partitions_file_format']['part_begin_col']]), int(row[input_config['partitions_file_format']['part_end_col']])))

	return partitions

def relevant_subpartitions(partitions, snp_first, snp_last):
	p_first = -1
	p_last = -1
	found_first = False
	found_last = False

	max_last = float('-inf')

	for i, p in enumerate(partitions):
		if p[0]<=snp_first and p[1]>=snp_first and not found_first:
			p_first = i
			found_first = True

		if p[0]<=snp_last and p[1]>=snp_last:
			p_last = i
			found_last = True

	if (not found_last) or (not found_first):
		raise Exception('Error: partition for snp_last or snp_first not found!')

	return partitions[p_first:(p_last+1)]

def first_last(name, input_config, first=-1, last=-1):
	'''
	evaluates first and last, and provides defaults if they are -1
	'''
	if first == -1 or last == -1:
		tmp_partitions = read_partitions(name, input_config)
		if first == -1:
			first = tmp_partitions[0][0]
		if last == -1:
			last = tmp_partitions[len(tmp_partitions)-1][1]

	return first, last

def insert_into_matrix(row, matrix, locus_list, partition_file_format, symmetric):
	loc_i = int(row[partition_file_format['i_col']])
	loc_j = int(row[partition_file_format['j_col']])

	if min((loc_i, loc_j)) == loc_i: 
		l = loc_i
		r = loc_j
		l_id = row[partition_file_format['i_id_col']] if partition_file_format['i_id_col'] >=0 else 'None'
		r_id = row[partition_file_format['j_id_col']] if partition_file_format['j_id_col'] >=0 else 'None'
		l_g = float(row[partition_file_format['g_i_col']])
		r_g = float(row[partition_file_format['g_j_col']])
	else:
		r = loc_i
		l = loc_j
		r_id = row[partition_file_format['i_id_col']] if partition_file_format['i_id_col'] >=0 else 'None'
		l_id = row[partition_file_format['j_id_col']] if partition_file_format['j_id_col'] >=0 else 'None'
		r_g = float(row[partition_file_format['g_i_col']])
		l_g = float(row[partition_file_format['g_j_col']])						

	naive = float(row[partition_file_format['naive_val_col']])
	shrink = float(row[partition_file_format['shrink_val_col']])

	if l not in matrix:
		matrix[l] = {}
		matrix[l]['data'] = {}
		matrix[l]['desc'] = {}
		matrix[l]['desc']['l_id'] = l_id
		matrix[l]['desc']['l_g'] = l_g
		# The SNPs in the input files are not guaranteed to be in sorted order!!
		where = bisect.bisect_left(locus_list, l)
		locus_list.insert(where, l)
		# locus_list.append(l)
		# heapq.heappush(locus_list, l)

	if r not in matrix[l]['data']:
		matrix[l]['data'][r] = {}
		matrix[l]['data'][r]['r_id'] = r_id
		matrix[l]['data'][r]['r_g'] = r_g
		matrix[l]['data'][r]['naive'] = naive
		matrix[l]['data'][r]['shrink'] = shrink
		# matrix[l]['data'][r]['corr_coeff'] = float('nan')

	if symmetric:
		if r not in matrix:
			matrix[r] = {}
			matrix[r]['data'] = {}
			matrix[r]['desc'] = {}
			matrix[r]['desc']['l_id'] = r_id
			matrix[r]['desc']['l_g'] = r_g
			# The SNPs in the input files are not guaranteed to be in sorted order!!
			where = bisect.bisect_left(locus_list, r)
			locus_list.insert(where, r)
			# locus_list.append(l)
			# heapq.heappush(locus_list, l)

		if l not in matrix[r]['data']:
			matrix[r]['data'][l] = {}
			matrix[r]['data'][l]['r_id'] = l_id
			matrix[r]['data'][l]['r_g'] = l_g
			matrix[r]['data'][l]['naive'] = naive
			matrix[r]['data'][l]['shrink'] = shrink

def insert_into_matrix_lean(row, matrix, locus_list, partition_file_format, symmetric):
	loc_i = int(row[partition_file_format['i_col']])
	loc_j = int(row[partition_file_format['j_col']])

	if min((loc_i, loc_j)) == loc_i: 
		l = loc_i
		r = loc_j
		l_id = row[partition_file_format['i_id_col']] if partition_file_format['i_id_col'] >=0 else 'None'
		r_id = row[partition_file_format['j_id_col']] if partition_file_format['j_id_col'] >=0 else 'None'
		l_g = float(row[partition_file_format['g_i_col']])
		r_g = float(row[partition_file_format['g_j_col']])
	else:
		r = loc_i
		l = loc_j
		r_id = row[partition_file_format['i_id_col']] if partition_file_format['i_id_col'] >=0 else 'None'
		l_id = row[partition_file_format['j_id_col']] if partition_file_format['j_id_col'] >=0 else 'None'
		r_g = float(row[partition_file_format['g_i_col']])
		l_g = float(row[partition_file_format['g_j_col']])						

	naive = float(row[partition_file_format['naive_val_col']])
	shrink = float(row[partition_file_format['shrink_val_col']])

	if l not in matrix:
		matrix[l] = {}
		# The SNPs in the input files are not guaranteed to be in sorted order!!
		where = bisect.bisect_left(locus_list, l)
		locus_list.insert(where, l)
		# locus_list.append(l)
		# heapq.heappush(locus_list, l)

	if r not in matrix[l]:
		matrix[l][r] = shrink
		# matrix[l]['data'][r]['corr_coeff'] = float('nan')

def read_partition_into_matrix(partitions, p_index, matrix, locus_list, name, input_config, snp_first, snp_last, symmetric=False):
	with gzip.open(input_config['partition_root']
					+name+'/'
					+name+input_config['partition_filename_delim']
					+str(partitions[p_index][0])+input_config['partition_filename_delim']
					+str(partitions[p_index][1])+input_config['partition_filename_ext'], 'rt') as f_in:
		csv_reader = csv.reader(f_in, delimiter=input_config['partition_delim'])
		try:
			for row in csv_reader:
					insert_into_matrix(row, matrix, locus_list, input_config['partition_file_format'], symmetric)

		except EOFError as e:
			print_log_msg('Error: '+str(e)+' at: '+str(partitions[p_index][0])+' '+str(partitions[p_index][1]))

def read_partition_into_matrix_lean(partitions, p_index, matrix, locus_list, name, input_config, snp_first, snp_last, symmetric=False):
	with gzip.open(input_config['partition_root']
					+name+'/'
					+name+input_config['partition_filename_delim']
					+str(partitions[p_index][0])+input_config['partition_filename_delim']
					+str(partitions[p_index][1])+input_config['partition_filename_ext'], 'rt') as f_in:
		csv_reader = csv.reader(f_in, delimiter=input_config['partition_delim'])
		try:
			for row in csv_reader:
					insert_into_matrix_lean(row, matrix, locus_list, input_config['partition_file_format'], symmetric)

		except EOFError as e:
			print_log_msg('Error: '+str(e)+' at: '+str(partitions[p_index][0])+' '+str(partitions[p_index][1]))

def write_output(out_filename, locus_list, locus_list_deleted, sum_list, out_delim, sum_list_len=None):
	with gzip.open(out_filename, 'wt') as f_out:
		csv_writer = csv.writer(f_out, delimiter=out_delim)

		print_log_msg('First part of output')

		print_log_msg('Size of locus_list_deleted: '+str(len(locus_list_deleted)))
		for l in locus_list_deleted:
			if l not in sum_list:
				print_log_msg(str(l)+' not in sum list!')				
			else:
				if sum_list_len is None:
					csv_writer.writerow([l, sum_list[l]])
				else:
					csv_writer.writerow([l, sum_list[l]/sum_list_len[l]])

		print_log_msg('Second part of output')

		print_log_msg('Size of locus_list: '+str(len(locus_list)))
		for l in locus_list:
			if l not in sum_list:
				print_log_msg(str(l)+' not in sum list!')				
			else:
				if sum_list_len is None:
					csv_writer.writerow([l, sum_list[l]])
				else:
					csv_writer.writerow([l, sum_list[l]/sum_list_len[l]])

	print_log_msg('Output done')

def write_output_lean(out_filename, locus_list, start, end, sum_list, out_delim, sum_list_len=None):
	with gzip.open(out_filename, 'at') as f_out:
		csv_writer = csv.writer(f_out, delimiter=out_delim)

		print_log_msg('Only part of output')

		print_log_msg('Size of locus_list: '+str(len(locus_list)))
		for i in range(start, end):
			if locus_list[i] not in sum_list:
				print_log_msg(str(locus_list[i])+' not in sum list!')				
			else:
				if sum_list_len is None:
					csv_writer.writerow([locus_list[i], sum_list[locus_list[i]]])
				else:
					csv_writer.writerow([locus_list[i], sum_list[locus_list[i]]/sum_list_len[locus_list[i]]])


	print_log_msg('Output done')

# def delete_loci_smaller_than_and_output_matrix_to_file(cutoff, matrix, locus_list, locus_list_deleted, out_fname, out_matrix_delim):
# 	print_log_msg('Writing matrix segment to output file')
# 	with gzip.open(out_fname, 'at') as f_out:
# 		csv_writer = csv.writer(f_out, delimiter=out_matrix_delim)
# 		cnt = 0
# 		while locus_list[cnt] < cutoff: # and cnt<len(locus_list)
# 			for key, val in matrix[locus_list[cnt]]['data'].items():
# 				if 'corr_coeff' in val:
# 					csv_writer.writerow( [locus_list[cnt], key, val['corr_coeff']] )
# 			cnt+=1

# 	delete_loci_smaller_than(cutoff, matrix, locus_list, locus_list_deleted)	

def delete_loci_smaller_than(cutoff, matrix, locus_list, locus_list_deleted):
	print_log_msg('Deleteing segment of matrix') 
	
	cnt = delete_matrix(cutoff, locus_list, matrix)
	
	locus_list_deleted.extend(locus_list[0:cnt])
	locus_list[0:cnt] = [] 

def delete_loci_smaller_than_leanest(cutoff, matrix, locus_list):
	print_log_msg('Deleteing segment of matrix, locus_list_deleted is not extended!') 
	
	cnt = delete_matrix(cutoff, locus_list, matrix)
	
	locus_list[0:cnt] = [] 

def delete_loci_smaller_than_lean(cutoff, matrix, locus_list, locus_list_deleted, out_fname, sum_list, out_delim):
	print_log_msg('Writing segment of matrix to file and deleting') 

	cnt = delete_matrix(cutoff, locus_list, matrix)
	# matrix = matrix.copy()
	
	# locus_list_deleted.extend(locus_list[0:cnt]) 
	write_output_lean(out_fname, locus_list, 0, cnt, sum_list, out_delim)
	for i in range(0, cnt):
		if locus_list[i] in sum_list:
			del sum_list[locus_list[i]]
	# sum_list = sum_list.copy()
	
	locus_list[0:cnt] = [] 

def delete_matrix(cutoff, locus_list, matrix):
	cnt = 0
	while cnt<len(locus_list) and locus_list[cnt] < cutoff: # and cnt<len(locus_list)
		del matrix[locus_list[cnt]]
		cnt+=1

	return cnt

# def determine_start_locus_vert(p_num, locus_list, snp_first):
# 	# Determine first locus
# 	curr_locus = -1
# 	if p_num == 0:
# 		if len(locus_list)>0:
# 			# Find first locus >= snp_first
# 			for i, locus in enumerate(locus_list):
# 				if locus >= snp_first:
# 					curr_locus = locus
# 					start_locus = locus
# 					curr_locus_index = i
# 					start_locus_index = i
# 					break
# 		else:
# 			raise Exception('Error: locus_list seems to be empty') 
# 	else:
# 		if len(locus_list)>0:
# 			curr_locus = locus_list[0]
# 			curr_locus_index = 0
# 		else:
# 			raise Exception('Error: locus_list seems to be empty')
# 
# 	if curr_locus<0:
# 		raise Exception('Error: curr_locus not found!')	
# 	
# 	return curr_locus, curr_locus_index, start_locus, start_locus_index
# 
# def determine_end_locus_vert(p_num, partitions, locus_list, snp_last):
# 	if p_num+1 < len(partitions):
# 		end_locus = partitions[p_num+1][0]
# 		end_locus_index = -1
# 	else:
# 		# end_locus = partitions[p_num][1]
# 	
# 		# Find last locus <= snp_last
# 		for i in reversed(range(0, len(locus_list))):
# 		# for locus in reversed(locus_list):
# 			if locus_list[i] <= snp_last:
# 				end_locus = locus_list[i]
# 				end_locus_index = i
# 				break
# 			
# 	return end_locus, end_locus_index
	

def print_log_msg(msg):
	print('['+time.strftime("%H:%M:%S")+'] '+msg) 
	sys.stdout.flush()

def read_hotspots(input_path_and_fname):
	with gzip.open(input_path_and_fname, 'rt') as f_in:
		graph_x = []
		graph_y = []
		pairs = []
		
		csv_reader = csv.reader(f_in, delimiter='\t')

		last_row = csv_reader.__next__() # skip 1st line
		last_row = csv_reader.__next__()
		graph_x.append(int(last_row[1]))
		graph_y.append(float(last_row[2]))

		for row in csv_reader:
			graph_x.append(int(row[1]))
			graph_y.append(float(last_row[2]))
			graph_x.append(int(row[1]))
			graph_y.append(float(row[2]))
			pairs.append((int(row[1]), float(row[2])))

			last_row = row

		return (graph_x,graph_y,pairs)

