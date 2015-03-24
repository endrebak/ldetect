import copy

tmp_path = '/nethome/jkpickrell/1kG_data/covariance_matrix/'

const = {}
const['orig_data'] = {
		'partition_root' : tmp_path,
		'partitions_dir' : tmp_path + 'scripts/',
		'partitions_file_suffix' : '_partitions',
		'partition_filename_delim' : '.',
		'partition_filename_ext' : '.gz',
		'partitions_delim' : ' ',
		'partition_delim' : ' ',
		'partition_file_format' : {
				'i_id_col':0, 
				'j_id_col':1, 
				'i_col':2, 
				'j_col':3, 
				'g_i_col':4, 
				'g_j_col':5, 
				'naive_val_col':6, 
				'shrink_val_col':7
		},
		'partitions_file_format' : {
				'part_begin_col':0,
				'part_end_col':1
		}
}

const['orig_data_EUR'] = copy.deepcopy(const['orig_data'])

tmp_path_ASN = '/nethome/jkpickrell/1kG_data/covariance_matrix/ASN/'
const['orig_data_ASN'] = copy.deepcopy(const['orig_data'])
const['orig_data_ASN']['partition_root'] =  tmp_path_ASN
const['orig_data_ASN']['partitions_dir'] = tmp_path_ASN + 'scripts/'

tmp_path_AFR = '/nethome/jkpickrell/1kG_data/covariance_matrix/AFR/'
const['orig_data_AFR'] = copy.deepcopy(const['orig_data'])
const['orig_data_AFR']['partition_root'] =  tmp_path_AFR
const['orig_data_AFR']['partitions_dir'] = tmp_path_AFR + 'scripts/' 

# const['part_buffer'] = 1
const['bin_size'] = 1

const['out_filename'] = 'output.txt.gz'
const['out_delim'] = '\t'

const['out_matrix_filename'] = 'output_matrix.txt.gz'
const['out_matrix_delim'] = '\t'

const['img_out_fname'] = 'output'
const['img_out_ext'] = '.png'

const['synth_input_root'] = '/data/research/pickrell_lab/genome_data/LD_covariance_matrix/ld_data/raw_synth_data/'
const['synth_input_filename_single'] = 'chrTest/ms_output.txt'


const['synth_data'] = copy.deepcopy(const['orig_data'])
const['synth_data']['partition_root'] = '/data/research/pickrell_lab/genome_data/LD_covariance_matrix/ld_data/covariance_matrix/'
const['synth_data']['partitions_dir'] = const['synth_data']['partition_root'] + 'scripts/'

tmp_path_maps = '/nethome/jkpickrell/Databases/human_genome/hg19/genetic_map/HapMap2_lifted/'

const['genetic_maps'] = {}
const['genetic_maps']['root'] = tmp_path_maps
const['genetic_maps']['file_base'] = 'genetic_map_GRCh37_'
const['genetic_maps']['ext'] = '.txt.gz'

const['gwas_plots'] = {}
const['gwas_plots']['height'] = {}
const['gwas_plots']['height']['root'] = '/data/research/pickrell_lab/genome_data/LD_covariance_matrix/height_gwas/'
const['gwas_plots']['height']['file_prefix'] = 'top00_noext_'
const['gwas_plots']['height']['file_suffix'] = '.bfs.bed'

