# notes

number entries all gz together: 244924474 (0.24e9)
number unique entries/snps: 372565
number snps in covar file: 196461

## todo

- ask about how to find recombination rates for missing populations

## possible bugs in original

- bug in calc_covariance? downstream scripts expect a position against itself to have a covariance, but it does not always have one

This is because a snp might be homogenous, i.e. all genotypes are equal in the
vcf for that SNP. This makes D zero, which means it is not included.

By using theta2 as a default value in matrix_to_vector this is fixed.

## possible bugs in mine

script that creates paritions makes multiple intervals with the same ends

18097260 18995988
18252321 18995988
18387002 18995988
18525697 18995988

## various

* original paper might have either used recombination rates from hapmap or CEU. Our method is better?

about 1000 people in 1kg sample info are not in the vcfs - what gives?

* in multiallelic sites, the alternate allele is not used for computation
