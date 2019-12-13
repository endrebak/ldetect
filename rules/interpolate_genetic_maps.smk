from glob import glob
import pandas as pd

chromosomes = ["chr" + str(i) for i in range(1, 23)]

prefix = config["prefix"]

f = "{prefix}/genetic_maps/{population}/{chromosome}.bed"
f = "{prefix}/genetic_maps/{population}/{chromosome}_hg38.bed"
f = "{prefix}/1kg/{chromosome}.vcf.gz.tbi"
f = "{prefix}/1kg/{chromosome}.vcf.gz"
f = "{prefix}/genetic_maps/{population}/{chromosome}_hg38.bed"

variant_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.{chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

f = "{prefix}/1kg/{chromosome}.vcf.gz.tbi"
f = "{prefix}/1kg/{chromosome}.vcf.gz"
# f = "{prefix}/genetic_maps/{population}/{chromosome}.txt.gz"

rule all:
    input:
        expand(f, chromosome=chromosomes, prefix=prefix, population="CEU")


rule download_genetic_maps:
    output:
        "{prefix}/genetic_maps/{population}.tar"
    shell:
        "curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/{wildcards.population}_omni_recombination_20130507.tar > {output[0]}"


rule untar_genetic_maps:
    input:
        "{prefix}/genetic_maps/{population}.tar"
    output:
        expand("{{prefix}}/genetic_maps/{{population}}/{chromosome}.txt.gz", chromosome=chromosomes)
    run:
        tmpfolder = wildcards.prefix + "/tmp/"
        shell("mkdir -p {tmpfolder}")
        shell("tar -xvf {input[0]} -C {tmpfolder}")

        for f in glob(tmpfolder + f"{wildcards.population}/**.gz", ):
            chromosome = "chr" + f.split("/")[-1].split("-")[1]
            final_file = f"{wildcards.prefix}/genetic_maps/{wildcards.population}/{chromosome}.txt.gz"
            shell("cp {f} {final_file}")

        shell("rm -rf {tmpfolder}")


rule genetic_maps_to_bed:
    input:
        "{prefix}/genetic_maps/{population}/{chromosome}.txt.gz"
    output:
        "{prefix}/genetic_maps/{population}/{chromosome}_hg19.bed"
    run:
        df = pd.read_table(input[0], sep="\t", usecols=[0, 2, 3], header=0)
        df.columns = "Start Name Score".split()

        df.insert(1, "End", df.Start + 1)
        df.insert(0, "Chromosome", "chr" + wildcards.chromosome)
        df.insert(5, "Strand", ".")

        df.to_csv(output[0], sep="\t", index=False, header=False)


rule fetchchain:
    output:
        "{prefix}/chain/hg19_to_hg38.chain.gz"
    shell:
        "curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz > {output[0]}"


rule hg19_to_38:
    input:
        chain = "{prefix}/chain/hg19_to_hg38.chain.gz",
        bed = "{prefix}/genetic_maps/{population}/{chromosome}_hg19.bed"
    output:
        bed = "{prefix}/genetic_maps/{population}/{chromosome}_hg38.bed",
        unmapped = "{prefix}/genetic_maps/{population}/{chromosome}.unmapped"
    shell:
        "liftOver {input.bed} {input.chain} {output.bed} {output.unmapped} && touch {output.unmapped}"


rule fetch_variants:
    output:
        "{prefix}/1kg/{chromosome}.vcf.gz"
    run:
        url = variant_url.format(chromosome=wildcards.chromosome)
        shell("curl {url} > {output[0]}")


rule fetch_variants_index:
    output:
        "{prefix}/1kg/{chromosome}.vcf.gz.tbi"
    run:
        url = variant_url.format(chromosome=wildcards.chromosome) + ".tbi"
        shell("curl {url} > {output[0]}")


rule vcf_to_bed:
    input:
        "{prefix}/1kg/{chromosome}.vcf.gz"
    output:
        "{prefix}/1kg/{chromosome}.bed"
    run:
        tmpfile = f"{wildcards.chromosome}.deleteme"

        shell("zcat {input[0]} | cut -f -2 > {tmpfile}")

        df = pd.read_table(tmpfile)
        df.columns = "Chromosome Start".split()
        df.insert(df.shape[1], "End", df.Start + 1)
        df.loc[:, "Chromosome"] = "chr" + df.Chromosome

        shell("rm {tmpfile}")

        df.to_csv(output[0], sep="\t")


rule interpolate_genetic_maps:
    input:
