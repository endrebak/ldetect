calc_covariance.py

Script used to calculate the Wen and Stephens shrinkage estimator of the covariance matrix. Reads a VCF file from stdin.

Input:

1. Genetic map file
2. Ouput covariance file  [gzipped]
3. List of individuals to use in the calculation

NOTE: two parameters are hard-coded. These are the effective population size (set to 11418, which is appropriate for European populations), and the cutoff below which the covariance is not reported (set to 1e-7). These can be changed in lines 10-11 in the script.


Example usage:

>tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr2.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 2:39967768-40367768 | python calc_covariance.py exampledata/chr2.interpolated_genetic_map.gz testout.gz exampledata/eurinds

This will output the covariance matrix to testout.gz. The columns are:

[snpid 1] [snpid 2] [position 1] [position 2] [genetic position 1] [genetic position 2] [empirical covariance] [shrinkage covariance]
