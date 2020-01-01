from glob import glob
import pandas as pd

# variant_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.{chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
variant_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"


rule download_genetic_maps:
    output:
        "{prefix}/genetic_maps/{population}.tar"
    shell:
        (
        "axel -q "
        "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/{wildcards.population}_omni_recombination_20130507.tar "
        "-o {output[0]}"
        )


rule untar_genetic_maps:
    input:
        rules.download_genetic_maps.output # "{prefix}/genetic_maps/{population}.tar"
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
        df.insert(0, "Chromosome", wildcards.chromosome)
        df.insert(5, "Strand", ".")

        df.to_csv(output[0], sep="\t", index=False, header=False)


rule fetchchain:
    output:
        "{prefix}/chain/hg19_to_hg38.chain.gz"
    shell:
        "curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz > {output[0]}"


rule genetic_maps_hg19_to_hg38:
    input:
        chain = rules.fetchchain.output, # "{prefix}/chain/hg19_to_hg38.chain.gz",
        bed = rules.genetic_maps_to_bed.output # "{prefix}/genetic_maps/{population}/{chromosome}_hg19.bed"
    output:
        bed = "{prefix}/genetic_maps/{population}/{chromosome}_hg38.bed",
        unmapped = "{prefix}/genetic_maps/{population}/{chromosome}.unmapped"
    shell:
        "liftOver {input.bed} {input.chain} {output.bed} {output.unmapped} && touch {output.unmapped}"


rule fetch_variants:
    output:
        "{prefix}/1kg/{chromosome}.vcf.gz"
    resources:
        instances = 1
    run:
        url = variant_url.format(chromosome=wildcards.chromosome)
        shell("axel {url} -q -o {output[0]}")



rule index_variants:
    input:
        rules.fetch_variants.output[0]
    output:
        "{prefix}/1kg/{chromosome}.vcf.gz.tbi"
    shell:
        "tabix {input[0]}"


rule vcf_to_bed:
    input:
        rules.fetch_variants.output
    output:
        "{prefix}/1kg/{chromosome}.bed"
    run:
        import gzip

        rows_to_skip = 0

        with gzip.open(input[0]) as fh:
            for i, l in enumerate(fh, 0):
                rows_to_skip = i
                if not l.decode().startswith("#"):
                    break

        print("rows_to_skip", rows_to_skip)

        df = pd.read_table(input[0], header=None, skiprows=rows_to_skip, usecols=[0, 1], nrows=None)
        print(df.head())
        df.columns = "Chromosome Start".split()
        df.insert(df.shape[1], "End", df.Start + 1)
        df.loc[:, "Chromosome"] = "chr" + df.Chromosome.astype(str)

        df.to_csv(output[0], sep="\t", index=False, header=False)


rule interpolate_genetic_maps:
    input:
        bed = rules.vcf_to_bed.output[0],
        mapfile = rules.genetic_maps_hg19_to_hg38.output.bed
    output:
        "{prefix}/genetic_maps/{population}/interpolated_{chromosome}_hg38.bed"
    script:
        "../scripts/interpolate_maps.py"
