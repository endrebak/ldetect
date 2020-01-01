import numpy as np

rule effective_population_size:
    input:
        config["effective_population_sizes"]
    output:
        "{prefix}/partitions/{population}/eps.txt"
    run:
        df = pd.read_table(input[0], header=0)
        egs = df[df.Population.isin(populations)].EGS.mean()

        with open(output[0], "w+") as egs_h:
            egs_h.write(str(egs) + "\n")




rule subset_on_population:
    input:
        variants = rules.fetch_variants.output[0],
        samples = rules.individuals_in_reference_panel.output.samples,
        index = rules.index_variants.output[0]
    output:
        protected("{prefix}/1kg/{population}/{chromosome}.vcf.gz")
    shell:
        "bcftools view --threads 48 --force-samples -O z -S {input.samples} {input.variants} > {output[0]}"


rule index_population_vcf:
    input:
        variants = rules.subset_on_population.output[0]
    output:
        "{prefix}/1kg/{population}/{chromosome}.vcf.gz.tbi"
    shell:
        "tabix -f {input[0]}"

covariate_outfile = "{prefix}/partition_covariances/{population}/{chromosome}/{start}_{end}.vcf.gz"
_covariate_outfile = "{prefix}/partition_covariances/{population}/{chromosome}/{{start}}_{{end}}.vcf.gz"

rule covariance_matrix:
    input:
        genetic_maps = "{prefix}/genetic_maps/CEU/interpolated_{chromosome}_hg38.bed",
        individuals_to_use = rules.individuals_in_reference_panel.output.samples,
        effective_population_size = rules.effective_population_size.output[0],
        variants =  rules.subset_on_population.output[0],
        index = rules.index_population_vcf.output[0]
    output:
        "{prefix}/partition_covariances/{population}/{chromosome}/{start}_{end}.tsv.gz"
    run:
        # df = pd.read_table(input.intervals, header=None, usecols=[0, 1], sep=" ", nrows=None)
        # df = df.tail(int(len(df)/2)).tail(1)
        # df = df.head(3).tail(1)
        # print(df)
        w = wildcards


        effective_population_size = -1
        with open(input.effective_population_size) as eps:
            effective_population_size = float(eps.readline().strip())

        covariance_cutoff = config["covariance_cutoff"]

        prefix = w.prefix
        population = w.population
        chromosome = w.chromosome
        start = w.start
        end = w.end

        c = chromosome.replace("chr", "")

        # outfile = _covariate_outfile.format(**locals())
        # outfile_template = output[0] + "/" + "{start}_{end}.txt.gz"

        template = f"""tabix {input.variants} {c}:{start}-{end} | cut -f 2,3,10- | tr "|" "\\t" |
python scripts/calc_covar.py {input.genetic_maps} {input.individuals_to_use} {effective_population_size} {covariance_cutoff} | gzip > {output[0]}"""
        # print(template)
        shell(template)

        # from time import time
        # for i, (_, (start, end)) in enumerate(df.iterrows()):
        #     tstart = time()
        #     print(i, len(df), i/len(df))
        #     outfile = outfile_template.format(start=start, end=end)
        #     cmd = template.format(start=start, end=end, outfile=outfile)
        #     shell(cmd)
        #     tend = time()
        #     print("Took:", tend - tstart)

rule intervals_to_parquet:
    input:
        "{prefix}/partition_covariances/{population}/{chromosome}/{start}_{end}.tsv.gz"
    output:
        "{prefix}/partition_covariances/{population}/{chromosome}/{start}_{end}.pq"
    run:
        f = input[0]
        df = pd.read_csv(f, sep=" ", usecols=[2, 3, 7], names="i j val".split(),
                         dtype={"i": np.int32, "j": np.int32})
        df.to_parquet(output[0])


def get_intervals(w):

    intervals = checkpoints.partition_chromosomes.get(**w).output[0]
    # print(intervals)
    df = pd.read_csv(intervals, header=None, sep=" ")
    # df = df.head(4)
    # print(df)
    starts = df[0].tolist()
    ends = df[1].tolist()
    f = rules.intervals_to_parquet.output
    pops = [w.population] * len(df)
    cs = [w.chromosome] * len(df)
    prefixes = [w.prefix] * len(df)
    fs = expand(f, zip, prefix=prefixes, population=pops, chromosome=cs, start=starts, end=ends)
    # print(fs)
    return fs

rule calculate_theta2:
    input:
        rules.individuals_in_reference_panel.output.samples
    output:
        "{prefix}/thetas2/{population}/{chromosome}.txt"
    run:
        inds = pd.read_table(input[0], header=None, squeeze=True).to_list()

        nind_int = len(inds)
        s = 0

        for _i in range(1, 2*nind_int):
            s = s+ 1.0/float(_i)

        s = 1/s

        theta = s/(2.0*float(nind_int)+s)
        thetas2 = (theta/2.0)*(1-theta/2.0)

        with open(output[0], "w+") as o:
            o.write(str(thetas2) + "\n")

rule collect_covariances:
    input:
        get_intervals
    output:
        "{prefix}/collected_covariances/{population}/{chromosome}.gz"
    shell:
        "zcat {input} | gzip -9 > {output[0]}"


# def glob_intervals(w):

#     checkpoint_output = checkpoints.covariance_matrix.get(**w).output[0]
#     print(checkpoint_output)

#     return checkpoint_output

# def theta2(f):
#     val = open(f).readline().strip()
#     print("theta2", val)
#     return double(val)


rule matrix_to_vector:
    input:
        covariances = get_intervals,
        partitions = "{prefix}/partitions/{population}/{chromosome}.gz",
        theta2 = rules.calculate_theta2.output[0]
        # checkpoints.partition_chromosomes.output
    output:
        "{prefix}/partitions/covariance/{population}/{chromosome}.gz"
    benchmark:
        "{prefix}/partitions/covariance/{population}/{chromosome}.txt"
    shell:
        "python scripts/matrix_to_vector.py {input.partitions} {input.theta2} {input.covariances} > {output[0]}"



rule matrix_to_vector_much_ram:
    input:
        covariances = get_intervals,
        partitions = "{prefix}/partitions/{population}/{chromosome}.gz",
        theta2 = rules.calculate_theta2.output[0]
        # checkpoints.partition_chromosomes.output
    output:
        "{prefix}/partitions/covariance_much_ram/{population}/{chromosome}.txt"
    benchmark:
        "{prefix}/partitions/covariance_much_ram/{population}/{chromosome}.bmark"
    shell:
        "python scripts/matrix_to_vector_much_ram.py {input.partitions} {input.theta2} {input.covariances} # > {output[0]}"


