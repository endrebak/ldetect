#!/usr/bin/env python3

import ldetect.baselib.flat_file_consts as cnst
import ldetect.baselib.flat_file as flat
import ldetect.baselib.binary_search as binsrch

import sys
import os.path
import math
import bisect

# import resource
import numpy as np

# def binary_search(a, x, lo=0, hi=None):   # can't use a to specify default for hi
#     hi = hi if hi is not None else len(a) # hi defaults to len(a)   
#     pos = bisect_left(a,x,lo,hi)          # find insertion position
#     return (pos if pos != hi and a[pos] == x else -1) # don't walk off the end

class MatrixAnalysis:
	def __init__(self, name, input_config, snp_first=-1, snp_last=-1):
		self.matrix = {}
		self.locus_list = []
		self.vert_sum = {}
		self.vert_sum_len = {}
		self.locus_list_deleted = []
		self.name = name
		self.input_config = input_config

		self.snp_first, self.snp_last = flat.first_last(name, input_config, snp_first, snp_last)

		self.partitions = flat.get_final_partitions(self.input_config, self.name, self.snp_first, self.snp_last)
		self.dynamic_delete = False
		self.calculation_complete = False
		self.start_locus = -1
		self.start_locus_index = -1
		self.end_locus = -1
		self.end_locus_index = -1

	def add_corr_coeff(self, corr_coeff, locus):
		if locus not in self.vert_sum:
			self.vert_sum[locus] = corr_coeff **2
			self.vert_sum_len[locus] = 1
		else:
			self.vert_sum[locus] += corr_coeff **2
			self.vert_sum_len[locus] += 1

	def generate_img(self, img_full_path, marked_snp=None):
		import numpy as np
		import matplotlib as mpl
		# mpl.use('svg')
		mpl.use('Agg')
		import matplotlib.pyplot as pt
		mpl.rcParams.update({'font.size': 22})
# 		import svgwrite

		# if center is None:
		# 	first = self.snp_first
		# 	last = self.snp_last
		# else:
		# 	if distance_in_snps is None:
		# 		raise Exception('Error: center is defined, but distance_in_snps is not!')
		# 	else:
		# 		first, last = self.query_locus_list(center, distance_in_snps)

		if not self.calculation_complete:
			raise Exception('Error: Calculation has not been completed prior to image generation')

		if self.dynamic_delete:
			raise Exception('Error: The matrix was dynamically deleted - cannot generate full image!')

		if len(self.matrix)<=0:
			raise Exception('Error: The matrix is emmpty or erroneous')

		flat.print_log_msg('Image init')	
		# svg_document = svgwrite.Drawing(filename = cnst.const['svg_out_fname'],
		#                                 size = (self.end_locus_index-self.start_locus_index, self.end_locus_index-self.start_locus_index))
		# Draw background 
		# svg_document.add(svg_document.rect(insert=(0, 0), size=('100%', '100%'), rx=None, ry=None, fill='rgb(0,0,0)'))

		plot_mtrx_size = self.end_locus_index-self.start_locus_index+1
		plot_mtrx = [[0 for x in range(plot_mtrx_size)] for x in range(plot_mtrx_size)]

		flat.print_log_msg('Plot matrix size: '+str(plot_mtrx_size))

		flat.print_log_msg('Matrix size: '+str(len(self.matrix)))

		flat.print_log_msg('locus_list size: '+str(len(self.locus_list)))

		flat.print_log_msg('locus_list_deleted size: '+str(len(self.locus_list_deleted)))

		x_values = [0 for x in range(plot_mtrx_size)]

		flat.print_log_msg('Generating image data')
		for loc_i in self.matrix:
			if loc_i >= self.snp_first and loc_i <= self.snp_last:
				x_values[self.locus_list.index(loc_i)-self.start_locus_index] = loc_i

				for loc_j in self.matrix[loc_i]['data']:
					if loc_j >= self.snp_first and loc_j <= self.snp_last:
						# if len(svg_loci)<svg_length:
						if 'corr_coeff' in self.matrix[loc_i]['data'][loc_j]:
							# color = 255* ( 1- ( self.matrix[loc_i]['data'][loc_j]['corr_coeff'] ** 2 ) ) 
							try:
								plot_mtrx[self.locus_list.index(loc_i)-self.start_locus_index][self.locus_list.index(loc_j)-self.start_locus_index] = ((self.matrix[loc_i]['data'][loc_j]['corr_coeff']) **2)
							except IndexError:
								print(self.locus_list.index(loc_i)-self.start_locus_index)
								print(len(plot_mtrx))
								print(self.locus_list.index(loc_j)-self.start_locus_index)
								print(len(plot_mtrx[self.locus_list.index(loc_i)-self.start_locus_index]))
							# svg_document.add(svg_document.rect(insert = (self.locus_list.index(loc_i)-self.start_locus_index, self.locus_list.index(loc_j)-self.start_locus_index),
					  #                              size = ('1', '1'), 
					  #                              fill = 'rgb(255,'+str(int(color))+','+str(int(color))+')'))
							# svg_loci.add(curr_locus)
						else:
							flat.print_log_msg("No 'corr_coef' key at: "+str(loc_i)+' '+str(loc_j))
							# raise Exception('WTF')

		flat.print_log_msg('Writing image file...')

		fig = pt.gcf()
		dpi = fig.get_dpi()
		fig_size = fig.get_size_inches()


		# pt.pcolor(np.array(plot_mtrx), cmap='Reds', vmin=0, vmax=1)
		pt.pcolormesh(np.array(plot_mtrx), cmap='binary', vmin=0, vmax=1)

		pt.colorbar()
		# x_values = np.array(x_values) # needs to be numpy array for pcolormesh()
		# X, Y = np.meshgrid(x_values, x_values)
		# pt.pcolormesh(X, Y, np.array(plot_mtrx), cmap='Reds', vmin=0, vmax=1)

		if marked_snp is not None:
			bpoint_loc = x_values.index(marked_snp)

			pt.scatter((bpoint_loc),(bpoint_loc), marker='x', color='green')

			flat.print_log_msg('SNP: '+repr(marked_snp)+' @ index: '+repr(bpoint_loc)+' in graph')

		fig = pt.gcf()
		fig.set_size_inches((40,30))

		pt.xlabel('SNP #')
		pt.ylabel('SNP #')
		pt.title('Correlation coefficient squared matrix') 

		pt.savefig(img_full_path)

		# return x_values
		# svg_document.save()

# 	def read_partitions(self, input_config):
# 		if self.snp_first>self.snp_last:
# 			raise Exception('Error: snp_first is larger than snp_last!')
# 
# 		flat.print_log_msg('Reading partitions file')
# 		partitions = flat.read_partitions(self.name, input_config)
# 
# 		flat.print_log_msg('Getting relevant partitions')
# 		partitions = flat.relevant_subpartitions(partitions, self.snp_first, self.snp_last)
# 
# 		if len(partitions)<=0:
# 			raise Exception('Error: There are no relevant subpartitions')
# 
# 		return partitions

	def write_output_to_file(self, filename, out_delim, avg=False):
		if not self.calculation_complete:
			raise Exception('Error: Calculation has not been completed prior to output file generation')

		flat.print_log_msg('Writing output to file')
		if avg:
			flat.write_output(filename, self.locus_list, self.locus_list_deleted, self.vert_sum, out_delim, self.vert_sum_len)
		else:
			flat.write_output(filename, self.locus_list, self.locus_list_deleted, self.vert_sum, out_delim)

	# def query_locus_list(self, locus, distance_in_snps):
	# 	'''
	# 	Returns SNP locations that are distance_in_snps away from locus (in both directions)

	# 	If distance_in_snps runs out of the range [snp_first, snp_last] -> return the snp_first or snp_last, depending on where it runs out
	# 	'''

	# 	# ind=self.locus_list.index(locus)

	# 	if locus < self.snp_first or locus > self.snp_last:
	# 		raise Exception('Provided locus is outside of range [snp_first,snp_last]')

	# 	# In case an  arbitrary locus is given (it doesn't necessarily have to be in the range)
	# 	ind=binsrch.find_le_ind(self.locus_list, locus)
		
	# 	ind_first=ind-distance_in_snps
	# 	if ind_first<0:
	# 		ind_first=0

	# 	ind_last=ind+distance_in_snps
	# 	if ind_last>=len(self.locus_list):
	# 		ind_last=len(self.locus_list)-1

	# 	out_snp_first = self.locus_list[ind_first]
	# 	if out_snp_first < self.snp_first:
	# 		out_snp_first = self.snp_first

	# 	out_snp_last = self.locus_list[ind_last] 
	# 	if out_snp_last > self.snp_last:
	# 		out_snp_last = self.snp_last 

	# 	return out_snp_first, out_snp_last

	def calc_vert(self, dynamic_delete=True, sum_both_sides=True):
		# flat.print_log_msg('Removing existing matrix output file')
		# try:
		#     os.remove(cnst.const['out_matrix_delim'])
		# except OSError:
		#     pass

		raise Exception('calc_vert is deprecated - check code before running!')

		self.dynamic_delete = dynamic_delete

		flat.print_log_msg('Start')

		for p_num, p in enumerate(self.partitions):
			flat.print_log_msg('Reading partition: '+str(p))
			flat.read_partition_into_matrix(self.partitions, p_num, self.matrix, self.locus_list, self.name, self.input_config, self.snp_first, self.snp_last)

			# Determine first locus
			curr_locus = -1
			if p_num == 0:
				if len(self.locus_list)>0:
					# Find first locus >= snp_first
					for i, locus in enumerate(self.locus_list):
						if locus >= self.snp_first:
							curr_locus = locus
							start_locus = locus
							curr_locus_index = i
							start_locus_index = i
							break
				else:
					raise Exception('Error: locus_list seems to be empty') 
			else:
				if len(self.locus_list)>0:
					curr_locus = self.locus_list[0]
					curr_locus_index = 0
				else:
					raise Exception('Error: locus_list seems to be empty')
		
			if curr_locus<0:
				raise Exception('Error: curr_locus not found!')
			
			if p_num+1 < len(self.partitions):
				end_locus = self.partitions[p_num+1][0]
				end_locus_index = -1
			else:
				# end_locus = partitions[p_num][1]
			
				# Find last locus <= snp_last
				for i in reversed(range(0, len(self.locus_list))):
				# for locus in reversed(locus_list):
					if self.locus_list[i] <= self.snp_last:
						end_locus = self.locus_list[i]
						end_locus_index = i
						break
			
			flat.print_log_msg('Running for partition: '+str(p))
			# This will not include the very last SNP of the complete range, but that shouldn't be too important since the end of the range shouldn't be a defining location for LD
			while curr_locus < end_locus: 
				for key, el in self.matrix[curr_locus]['data'].items():
					corr_coeff = self.matrix[curr_locus]['data'][key]['shrink'] / math.sqrt( self.matrix[curr_locus]['data'][curr_locus]['shrink'] * self.matrix[key]['data'][key]['shrink'] )
					self.add_corr_coeff(corr_coeff, curr_locus)
					if sum_both_sides:
						self.add_corr_coeff(corr_coeff, key)

					# Just save it in the matrix ;)
					self.matrix[curr_locus]['data'][key]['corr_coeff'] = corr_coeff

				if curr_locus_index+1 < len(self.locus_list):
					curr_locus_index+=1
					curr_locus = self.locus_list[curr_locus_index]
				else:
					flat.print_log_msg('curr_locus_index out of bounds')
					break

			# flat.delete_loci_smaller_than_and_output_matrix_to_file(end_locus, self.matrix, locus_list, locus_list_deleted, cnst.const['out_matrix_filename'])
			if self.dynamic_delete:
				flat.print_log_msg('Deleting loci not required any more')
				flat.delete_loci_smaller_than(end_locus, self.matrix, self.locus_list, self.locus_list_deleted)

		self.start_locus = start_locus
		self.start_locus_index = start_locus_index
		self.end_locus = end_locus
		self.end_locus_index = end_locus_index

		self.calculation_complete = True

	def calc_diag_lean(self, out_fname, out_delim, dynamic_delete=True): 
		# flat.print_log_msg('Removing existing matrix output file')
		# try:
		#     os.remove(cnst.const['out_matrix_delim'])
		# except OSError:
		#     pass

		if dynamic_delete == False:
			raise Exception('Error: Conversion has been run in lean mode, but with dynamically=False.')

		self.dynamic_delete = dynamic_delete

		flat.print_log_msg('Start')

        # pre-read all relevant partitions at beginning!
		last_p_num = -1
		for p_num_init in range(0, len(self.partitions)-1):
			if self.snp_first >= self.partitions[p_num_init+1][0]:
				flat.print_log_msg('Pre-reading partition: '+str(self.partitions[p_num_init])) 
				flat.read_partition_into_matrix_lean(self.partitions, p_num_init, self.matrix, self.locus_list, self.name, self.input_config, self.snp_first, self.snp_last)
				last_p_num = p_num_init
			else:
				break

		curr_locus = -1
		# for p_num, p in enumerate(self.partitions):
		for p_num in range(last_p_num+1, len(self.partitions)):
			p = self.partitions[p_num]

			flat.print_log_msg('Reading partition: '+str(p))
			flat.read_partition_into_matrix_lean(self.partitions, p_num, self.matrix, self.locus_list, self.name, self.input_config, self.snp_first, self.snp_last)

			# Determine first locus
			if curr_locus<0: # Either first partition or not found in first partition
				# curr_locus = -1 # <- this should have been set to -1 before entering the main for loop
				if len(self.locus_list)>0:
					# Find first locus >= snp_first
					for i, locus in enumerate(self.locus_list):
						if locus >= self.snp_first:
							curr_locus = locus
							start_locus = locus
							curr_locus_index = i
							start_locus_index = i
							break
				else:
					raise Exception('Error: locus_list seems to be empty') 
			# else:
			# 	if len(self.locus_list)>0:
			# 		curr_locus = self.locus_list[0]
			# 		curr_locus_index = 0
			# 	else:
			# 		raise Exception('Error: locus_list seems to be empty')
			else:
				try:
					curr_locus_index = self.locus_list.index(curr_locus)
					# curr_locus is carried from prev iteration, but index has changed since part of matrix (and locus_list) has been deleted
				except ValueError:
					if len(self.locus_list)>0:
						curr_locus = self.locus_list[0]
						curr_locus_index = 0
					else:
						raise Exception('Error: locus_list seems to be empty')

			if curr_locus<0:
				flat.print_log_msg('Warning: curr_locus not found! Continuing to next partition.')
				flat.print_log_msg('Comment: This is possibly due to snp_first being very close to end of partition.')
				flat.print_log_msg('Details: ')
				flat.print_log_msg('Partition: '+repr(p))
				flat.print_log_msg('snp_first: '+repr(self.snp_first))
				flat.print_log_msg('curr_locus: '+repr(curr_locus)) 
				continue #continue to next partition 
				# raise Exception('Error: curr_locus not found!')	

			# Determine end locus
			if p_num+1 < len(self.partitions):
				end_locus = int((self.partitions[p_num][1] + self.partitions[p_num+1][0]) / 2) # diag - specific
			else:
				# end_locus = self.partitions[p_num][1]

				# Find last locus <= snp_last
				end_locus_found = False
				for i in reversed(range(0, len(self.locus_list))):
				# for locus in reversed(locus_list):
					if self.locus_list[i] <= self.snp_last:
						end_locus = self.locus_list[i]
						end_locus_index = i
						end_locus_found = True
						break

				if not end_locus_found:
					end_locus_index = 0
					end_locus = self.locus_list[end_locus_index]

			flat.print_log_msg('Running for partition: '+str(p))
			# This will not include the very last SNP of the complete range, but that shouldn't be too important since the end of the range shouldn't be a defining location for LD
			while curr_locus <= end_locus: 
				x = self.locus_list[curr_locus_index]
				y = self.locus_list[curr_locus_index]
				delta = 0

				while x >= self.partitions[p_num][0] and y <= self.partitions[p_num][1]:
					if x in self.matrix and y in self.matrix[x]:
						corr_coeff = self.matrix[x][y] / math.sqrt( self.matrix[x][x] * self.matrix[y][y] )
						self.add_corr_coeff(corr_coeff, curr_locus)
						# Just save it in the matrix ;) - removed for chrom11
						# self.matrix[x]['data'][y]['corr_coeff'] = corr_coeff
					# else:
					# 	flat.print_log_msg('Condition not satisfied 1!')
					# 	flat.print_log_msg('x: '+repr(x)+' y: '+repr(y))

					if delta!=0:
						x = self.locus_list[curr_locus_index-delta+1]
						if x in self.matrix and y in self.matrix[x]:
							corr_coeff = self.matrix[x][y] / math.sqrt( self.matrix[x][x] * self.matrix[y][y] )
							self.add_corr_coeff(corr_coeff, curr_locus)
							# Just save it in the matrix ;) - removed for chrom11
							# self.matrix[x]['data'][y]['corr_coeff'] = corr_coeff
						# else:
						# 	flat.print_log_msg('Condition not satisfied 2!')
						# 	flat.print_log_msg('x: '+repr(x)+' y: '+repr(y))

					delta += 1
					if curr_locus_index-delta >= 0:
						x = self.locus_list[curr_locus_index-delta]
					else:
						# flat.print_log_msg('X index out of bounds')
						break

					if curr_locus_index+delta < len(self.locus_list):
						y = self.locus_list[curr_locus_index+delta]
					else:
						# flat.print_log_msg('Y index out of bounds')
						break						

				if curr_locus_index+1 < len(self.locus_list):
					curr_locus_index+=1
					curr_locus = self.locus_list[curr_locus_index]
				else:
					flat.print_log_msg('curr_locus_index out of bounds')
					break

			# flat.print_log_msg('Mem before delete: '+repr(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
			# flat.delete_loci_smaller_than_and_output_matrix_to_file(end_locus, self.matrix, locus_list, locus_list_deleted, cnst.const['out_matrix_filename'])
			if self.dynamic_delete:
				flat.print_log_msg('Deleting loci not required any more')
				if p_num+1 < len(self.partitions):
					delete_loc = self.partitions[p_num+1][0]
				else:
					delete_loc = end_locus

				flat.delete_loci_smaller_than_lean(delete_loc, self.matrix, self.locus_list, self.locus_list_deleted, out_fname, self.vert_sum, out_delim)
			else:
				flat.print_log_msg('locus_list size: '+repr(len(self.locus_list)))

			# flat.print_log_msg('Mem after delete:  '+repr(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

		self.start_locus = start_locus
		self.start_locus_index = start_locus_index
		self.end_locus = end_locus
		self.end_locus_index = end_locus_index

		self.calculation_complete = True

	def calc_diag(self, dynamic_delete=True):
		# flat.print_log_msg('Removing existing matrix output file')
		# try:
		#     os.remove(cnst.const['out_matrix_delim'])
		# except OSError:
		#     pass

		self.dynamic_delete = dynamic_delete

		flat.print_log_msg('Start')

		# pre-read all relevant partitions at beginning!
		last_p_num = -1
		for p_num_init in range(0, len(self.partitions)-1):
			if self.snp_first >= self.partitions[p_num_init+1][0]:
				flat.print_log_msg('Pre-reading partition: '+str(self.partitions[p_num_init])) 
				flat.read_partition_into_matrix(self.partitions, p_num_init, self.matrix, self.locus_list, self.name, self.input_config, self.snp_first, self.snp_last)
				last_p_num = p_num_init
			else:
				break

		curr_locus = -1
		# for p_num, p in enumerate(self.partitions):
		for p_num in range(last_p_num+1, len(self.partitions)):
			p = self.partitions[p_num]

			flat.print_log_msg('Reading partition: '+str(p))
			flat.read_partition_into_matrix(self.partitions, p_num, self.matrix, self.locus_list, self.name, self.input_config, self.snp_first, self.snp_last)

			# Determine first locus
			if curr_locus<0: # Either first partition or not found in first partition
				# curr_locus = -1 # <- this should have been set to -1 before entering the main for loop
				if len(self.locus_list)>0:
					# Find first locus >= snp_first
					for i, locus in enumerate(self.locus_list):
						if locus >= self.snp_first:
							curr_locus = locus
							start_locus = locus
							curr_locus_index = i
							start_locus_index = i
							break
				else:
					raise Exception('Error: locus_list seems to be empty') 
			# else:
			# 	if len(self.locus_list)>0:
			# 		curr_locus = self.locus_list[0]
			# 		curr_locus_index = 0
			# 	else:
			# 		raise Exception('Error: locus_list seems to be empty')
			else:
				try:
					curr_locus_index = self.locus_list.index(curr_locus)
					# curr_locus is carried from prev iteration, but index has changed since part of matrix (and locus_list) has been deleted
				except ValueError:
					if len(self.locus_list)>0:
						curr_locus = self.locus_list[0]
						curr_locus_index = 0
					else:
						raise Exception('Error: locus_list seems to be empty')

			if curr_locus<0:
				flat.print_log_msg('Warning: curr_locus not found! Continuing to next partition.')
				flat.print_log_msg('Comment: This is possibly due to snp_first being very close to end of partition.')
				flat.print_log_msg('Details: ')
				flat.print_log_msg('Partition: '+repr(p))
				flat.print_log_msg('snp_first: '+repr(self.snp_first))
				flat.print_log_msg('curr_locus: '+repr(curr_locus)) 
				continue #continue to next partition
				# raise Exception('Error: curr_locus not found!')	
			

			# Determine end locus
			if p_num+1 < len(self.partitions):
				end_locus = int((self.partitions[p_num][1] + self.partitions[p_num+1][0]) / 2)
			else:
				# end_locus = self.partitions[p_num][1]

				# Find last locus <= snp_last
				end_locus_found = False
				for i in reversed(range(0, len(self.locus_list))):
				# for locus in reversed(locus_list):
					if self.locus_list[i] <= self.snp_last:
						end_locus = self.locus_list[i]
						end_locus_index = i
						end_locus_found = True
						break

				if not end_locus_found:
					end_locus_index = 0
					end_locus = self.locus_list[end_locus_index]

			flat.print_log_msg('Running for partition: '+str(p))
			# This will not include the very last SNP of the complete range, but that shouldn't be too important since the end of the range shouldn't be a defining location for LD
			while curr_locus <= end_locus: 
				x = self.locus_list[curr_locus_index]
				y = self.locus_list[curr_locus_index]
				delta = 0

				while x >= self.partitions[p_num][0] and y <= self.partitions[p_num][1]:
					if x in self.matrix and y in self.matrix[x]['data']:
						corr_coeff = self.matrix[x]['data'][y]['shrink'] / math.sqrt( self.matrix[x]['data'][x]['shrink'] * self.matrix[y]['data'][y]['shrink'] )
						self.add_corr_coeff(corr_coeff, curr_locus)

						# Just save it in the matrix ;) ...for img
						self.matrix[x]['data'][y]['corr_coeff'] = corr_coeff

					if delta!=0:
						x = self.locus_list[curr_locus_index-delta+1]
						if x in self.matrix and y in self.matrix[x]['data']:
							corr_coeff = self.matrix[x]['data'][y]['shrink'] / math.sqrt( self.matrix[x]['data'][x]['shrink'] * self.matrix[y]['data'][y]['shrink'] )
							self.add_corr_coeff(corr_coeff, curr_locus)

							# Just save it in the matrix ;) ...for img
							self.matrix[x]['data'][y]['corr_coeff'] = corr_coeff


					delta += 1
					if curr_locus_index-delta >= 0:
						x = self.locus_list[curr_locus_index-delta]
					else:
						# flat.print_log_msg('X index out of bounds')
						break

					if curr_locus_index+delta < len(self.locus_list):
						y = self.locus_list[curr_locus_index+delta]
					else:
						# flat.print_log_msg('Y index out of bounds')
						break						

				if curr_locus_index+1 < len(self.locus_list):
					curr_locus_index+=1
					curr_locus = self.locus_list[curr_locus_index]
				else:
					flat.print_log_msg('curr_locus_index out of bounds')
					break

			# flat.delete_loci_smaller_than_and_output_matrix_to_file(end_locus, self.matrix, locus_list, locus_list_deleted, cnst.const['out_matrix_filename'])
			if self.dynamic_delete:
				flat.print_log_msg('Deleting loci not required any more')
				if p_num+1 < len(self.partitions):
					delete_loc = self.partitions[p_num+1][0] # diag - specific
				else:
					delete_loc = end_locus

				flat.delete_loci_smaller_than(delete_loc, self.matrix, self.locus_list, self.locus_list_deleted)
			else:
				flat.print_log_msg('locus_list size: '+repr(len(self.locus_list)))

		self.start_locus = start_locus
		self.start_locus_index = start_locus_index
		self.end_locus = end_locus
		self.end_locus_index = end_locus_index

		self.calculation_complete = True 


# main() functionality has moved to the pipeline file: P02_matrix_to_vector_pipeline.py

# def main():
# 	if(sys.argv[6]=='img-yes'):
# 		generate_img = True
# 	elif(sys.argv[6]=='img-no'):
# 		generate_img = False
# 	else:
# 		raise Exception('Error: Unknown argument: '+sys.argv[6])

# 	analysis = MatrixAnalysis(sys.argv[1], cnst.const[sys.argv[4]], int(sys.argv[2]), int(sys.argv[3]))

# 	if(sys.argv[7]=='vert'):
# 		analysis.calc_vert(not generate_img) 
# 	elif(sys.argv[7]=='diag'):
# 		analysis.calc_diag(not generate_img) 
# 	else:
# 		raise Exception('Error: Unknown argument: '+sys.argv[7])

# 	if(sys.argv[5]=='avg'):
# 		avg = True
# 		raise Exception('Average used, but its output is not always consistent - especially for diag!')
# 	elif(sys.argv[5]=='sum'):
# 		avg = False
# 	else:
# 		raise Exception('Error: Unknown argument: '+sys.argv[5])

# 	analysis.write_output_to_file(cnst.const['out_filename'], cnst.const['out_delim'], avg)

# 	if generate_img:
# 		analysis.generate_img(cnst.const['img_out_fname']+cnst.const['img_out_ext'])

# 	flat.print_log_msg('Done')

# if __name__ == '__main__':
# 	main()
