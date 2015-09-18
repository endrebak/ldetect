# LDetect

## Input data

### Genetic maps

Genetic maps are available at:
https://github.com/joepickrell/1000-genomes-genetic-maps

### Reference panel

Reference panel example:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/[ALL.chr2.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz]

## Step 1 - Partition chromosome(s)

The idea is to split the chromosomes up into partitions by virtue of having large genetic distances (for example, so that they can run in parallel)

Input:

1. Genetic map
2. Number of individuals in reference panel

Usage:

python3 P00_00_partition_chromosome.py <input_genetic_map> <n_individuals_in_ref_panel> <output_file>

Example:

python3 P00_00_partition_chromosome.py example_data/chr2.interpolated_genetic_map.gz 379 example_data/cov_matrix/scripts/chr2_partitions

Output:

This will output the partitions to chr2_partitions. The columns are:
[start position] [stop position]

Caveat:

If the output directory for the covariance matrix is:
example_data/cov_matrix/
then downstream scripts EXPECT partition files to be in the scripts/ subdirectory and to be named chr<chr_number>_partitions

## Step 2 - Calculate covariance matrix

Script used to calculate the Wen and Stephens shrinkage estimator of the covariance matrix. Reads a VCF file from stdin.

Input:

1. Reference panel (via stdin)
2. Genetic map file
3. Ouput covariance file [gzipped]
4. List of individuals to use in the calculation
5. Effective population size
6. Cutoff beneath which covariance is not reported

Usage:

tabix -h <input_ref_panel> <chr_number>:<start>-<stop> | python3 P00_01_calc_covariance.py <input_genetic_map> <input_individuals_file> <effective_population_size> <cov_cutoff> <output_cov_matrix_partition>

Example:

tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr2.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 2:39967768-40067768 | python3 P00_01_calc_covariance.py example_data/chr2.interpolated_genetic_map.gz example_data/eurinds.txt 11418 1e-7 example_data/cov_matrix/chr2/chr2.39967768.40067768.gz

Output:

This will output the covariance matrix to chr2.39967768.40067768.gz. The columns are:

[snpid 1] [snpid 2] [position 1] [position 2] [genetic position 1] [genetic position 2] [empirical covariance] [shrinkage covariance]

NOTE: The effective population size in the example is set to 11418, which is appropriate for European populations, and the cutoff below which the covariance is not reported was set to 1e-7 in previous work.

Caveat:

If the output directory for the covariance matrix is:
example_data/cov_matrix/
then downstream scripts EXPECT covariance matrix files to be in subdirectories named chr<chr_number>/ and to be named chr<chr_number>.<start>.<stop>.gz, where <start> and <stop> correspond to values from the partitions file
