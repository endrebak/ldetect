import ldetect.baselib.flat_file_consts as cnst
import ldetect.baselib.flat_file as flat

import gzip
import os
import random
import subprocess

def ms_to_flat(name, input_path_and_fname, output_config, init_loc):
	with open(input_path_and_fname, 'r') as f_in:
		ignore = f_in.readline()
		ignore = f_in.readline()
		ignore = f_in.readline()

		curr_loc = init_loc

		while str.strip(f_in.readline()) == '//':
			segsites = f_in.readline()
			numsites = int(segsites.split(' ')[1])
			ignore = f_in.readline()
			ignore = f_in.readline()

			lines = []

			line = str.strip(f_in.readline())
			averages = [0] * len(line)

			# read from file
			while line != '':
				for i, c in enumerate(line):
					averages[i] += float(line[i])
				lines.append(line)
				line = str.strip(f_in.readline())

			for i in range(0, len(averages)):
				averages[i] /= len(lines)

			submatrix = []

			for i in range(0,len(lines[0])):
				col = []
				for j in range(i, len(lines[0])):
					val = 0
					for line in lines:
						val += (float(line[i])-averages[i])*(float(line[j])-averages[j])
						
					val /= (len(lines)-1)
					col.append(val)

				submatrix.append(col)

			end_loc = curr_loc+len(lines[0])-1

			append_partitions_file(name, curr_loc, end_loc, output_config)

			create_partition_file(name, curr_loc, end_loc, submatrix, output_config)

			curr_loc = end_loc+1

	return curr_loc

def append_partitions_file(name, begin, end, output_config):
	with open(output_config['partitions_dir']+name+output_config['partitions_file_suffix'], 'at') as f_out:
		f_out.write(str(begin)+output_config['partitions_delim']+str(end)+'\n')

def create_partition_file(name, begin, end, submatrix, output_config):
	with gzip.open(
				output_config['partition_root']+name+'/'+
				name+output_config['partition_filename_delim']+
				str(begin)+output_config['partition_filename_delim']+
				str(end)+output_config['partition_filename_ext']
			, 'wt') as f_out:
		for i in range(0,len(submatrix)):
			for j in range(0, len(submatrix[i])):
				loc1 = begin+i
				loc2 = begin+i+j
				name1 = 'rs-'+str(loc1)
				name2 = 'rs-'+str(loc2)

				f_out.write(	
					name1+output_config['partition_delim']+
					name2+output_config['partition_delim']+
					str(loc1)+output_config['partition_delim']+
					str(loc2)+output_config['partition_delim']+
					'0.0'+output_config['partition_delim']+
					'0.0'+output_config['partition_delim']+
					'0.0'+output_config['partition_delim']+
					str(submatrix[i][j])+'\n'
				)

def multiple_ms_to_flat(name, file_list, output_config, init_loc):
	os.mkdir(output_config['partition_root']+name)

	partitions_fname = output_config['partitions_dir']+name+output_config['partitions_file_suffix']
	if os.path.isfile(partitions_fname):
		raise Exception("Partitions file already exists!\n"+partitions_fname)

	loc = 1
	for fname in file_list:
		loc = ms_to_flat(name, fname, cnst.const['synth_data'], loc)

def generate_random_hotspots(all_hotspots, n_sites):
	if n_sites <= 0:
		raise Exception('Error: n_sites <= 0')

	last_hotspot_loc = all_hotspots[len(all_hotspots)-1][0]

	last_possible_index = 0
	for i in reversed(range(0,len(all_hotspots))):
		if last_hotspot_loc - all_hotspots[i][0] > n_sites:
			last_possible_index = i
			break

	first_index = random.randint(0,last_possible_index)

	chosen_hotspots = []
	for i in range(first_index, len(all_hotspots)):
		chosen_hotspots.append(all_hotspots[i])
		if all_hotspots[i][0] - all_hotspots[first_index][0] > n_sites:
			break

	return chosen_hotspots

def create_cmd_line_hotspots(chosen_hotspots, n_sites): # -> calculate true average! -> then relative -> then regions
	truncated_hotspots = chosen_hotspots
	truncated_hotspots[len(truncated_hotspots)-1] = (truncated_hotspots[0][0]+n_sites, float('nan')) # this just gives us a placeholder for the last boundary

	weighted_sum = 0
	for i in range(1, len(truncated_hotspots)):
		weighted_sum += (truncated_hotspots[i][0]-truncated_hotspots[i-1][0])*truncated_hotspots[i-1][1]

	wieghted_avg = weighted_sum / n_sites

	cmd_line_hotspots = ' '
	shift = truncated_hotspots[0][0]-1 # The hotspots should start from 1, not 0
	num_args = 0
	for i in range(1, len(truncated_hotspots)):
		begin = truncated_hotspots[i-1][0]-shift
		end = truncated_hotspots[i][0]-shift-1
		if begin!=end: # msHOT does not like regions starting and ending at the same locus
			cmd_line_hotspots += str(begin)+' '+str(end)+' '+str(truncated_hotspots[i-1][1]/wieghted_avg)+' '
			num_args += 1

	cmd_line_hotspots = str(num_args)+' '+cmd_line_hotspots # number of hotspots is first argument
	return cmd_line_hotspots

def run_multiple_msHOT(in_name, genetic_maps, out_name, output_root, num_runs, msHOT_params):
	out_dir_name = output_root+out_name+'/'
	os.mkdir(out_dir_name)

	#read in hotspot map
	trash_x, trash_y, all_hotspots = flat.read_hotspots(genetic_maps['root']+genetic_maps['file_base']+in_name+genetic_maps['ext'])

	# chosen_hotspots = generate_random_hotspots(all_hotspots, msHOT_params['total_n_sites'])
	# print(chosen_hotspots[len(chosen_hotspots)-1][0]-chosen_hotspots[0][0])

	file_list = []
	for file_num in range(0, num_runs):
		chosen_hotspots = generate_random_hotspots(all_hotspots, msHOT_params['total_n_sites'])
		cmd_line_hotspots = create_cmd_line_hotspots(chosen_hotspots, msHOT_params['total_n_sites']) # -> calculate true average! -> then relative -> then regions
		curr_file_name = out_dir_name+out_name+'-'+str(file_num)
		# run msHOT w/ params > curr_file_name
		subprocess.call('msHOT '+str(msHOT_params['sample_size'])+' '+str(msHOT_params['n_runs'])+' -t '+str(msHOT_params['theta'])+' -r '+str(msHOT_params['bgrd_rho'])+' '+str(msHOT_params['total_n_sites'])+' -v '+cmd_line_hotspots+' > '+curr_file_name, shell=True) 
		file_list.append(curr_file_name)

	return file_list

def main():
	input_chr_name = 'chr1'
	# output_chr_name = 'chrTestHOT_1run_100k'
	# num_runs = 1
	num_runs = 20

	N = 1e4 # population size
	mu = 1.5e-8 # mutation rate [/site]
	r = 1.5e-8 # recombination rate [/site] 
	# chrom_length = 1e5 # total chrom length
	chrom_length = 2e6

	# ms 400 1 -t 60 -r 60 100000 -v ...
	msHOT_params = {
		'sample_size' : 400,
		'n_runs' : 1,
		'theta' : 4*N*mu*chrom_length,
		'bgrd_rho' : 4*N*r*chrom_length,
		'total_n_sites' : chrom_length
	}

	output_chr_name = 'chrTestHOT_'+str(num_runs)+'runs_'+str(msHOT_params['total_n_sites']/1e6)+'M'

	msHOT_file_list = run_multiple_msHOT(input_chr_name, cnst.const['genetic_maps'], output_chr_name, cnst.const['synth_input_root'], num_runs, msHOT_params)

	# multiple_ms_to_flat(output_chr_name, cnst.const['synth_input_root']+cnst.const['synth_input_filename_single'], cnst.const['synth_data'], 1)
	multiple_ms_to_flat(output_chr_name, msHOT_file_list, cnst.const['synth_data'], 1)

if __name__ == '__main__':
	main()