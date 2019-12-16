import pandas as pd

sample_sheet_path = config["sample_sheet"]
prefix = config["prefix"]
print(sample_sheet_path)
ss = pd.read_table(sample_sheet_path)

print(ss)

rule all:
    input:
        expand("{prefix}/partitions/{chromosome}.gz", prefix=prefix, chromosome=ss.Chromosome)


rule download_chromosomes_vcf:
    output:
        "{prefix}/download/{chromosome}.gz"
    params:
        url = lambda w: ss[ss.Chromosome == w.chromosome].URL.iloc[0]
    shell:
        "curl {params.url} > {output[0]}"


rule liftover_recombination_maps:
    input:
        config["genetic_map"]
    output:
        "{prefix}/download/{chromosome}.gz"


rule partition_chromosomes:
    input:
        genetic_map = config["genetic_map"],
        number_individuals_in_reference_panel = ""
    output:
        "{prefix}/partitions/{chromosome}.gz"
    # params:
        
    script:
        "scripts/partition_chromosome.py"


