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


covariate_outfile = "{prefix}/partitions/{population}/{chromosome}/{start}_{end}.vcf.gz"
_covariate_outfile = "{prefix}/partitions/{population}/{chromosome}/{{start}}_{{end}}.vcf.gz"


rule subset_on_population:
    input:
        variants = rules.fetch_variants.output[0],
        samples = rules.individuals_in_reference_panel.output.samples,
        index = rules.fetch_variants_index.output[0]
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



checkpoint covariance_matrix:
    input:
        genetic_maps = "{prefix}/genetic_maps/CEU/interpolated_{chromosome}_hg38.bed",
        individuals_to_use = rules.individuals_in_reference_panel.output.samples,
        effective_population_size = rules.effective_population_size.output[0],
        variants =  rules.subset_on_population.output[0],
        intervals = rules.partition_chromosomes.output[0],
        index = rules.index_population_vcf.output[0]
    output:
        directory("{prefix}/partition_covariances/{population}/{chromosome}/")
    run:
        df = pd.read_table(input.intervals, header=None, usecols=[0, 1], sep=" ", nrows=None)
        df = df.tail(int(len(df)/2)).tail(1)
        # df = df.head(3).tail(1)
        print(df)
        w = wildcards

        effective_population_size = -1
        with open(input.effective_population_size) as eps:
            effective_population_size = float(eps.readline().strip())

        covariance_cutoff = config["covariance_cutoff"]

        prefix = w.prefix
        population = w.population
        chromosome = w.chromosome

        c = chromosome.replace("chr", "")

        outfile = _covariate_outfile.format(**locals())

        template = f"""tabix {input.variants} {c}:{{start}}-{{end}} | cut -f 2,3,10- | tr "|" "\\t" |
python scripts/calc_covar.py {input.genetic_maps} {input.individuals_to_use} {effective_population_size} {covariance_cutoff} {outfile}"""

        print(df.shape)
        for i, (_, (start, end)) in enumerate(df.iterrows()):
            print(i)
            cmd = template.format(start=start, end=end)
            print(cmd)
            shell(cmd)


def glob_intervals(w):

    checkpoint_output = checkpoints.covariance_matrix.get(**w).output[0]
    print(checkpoint_output)

    return checkpoint_output


rule matrix_to_vector:
    input:
        glob_intervals
    output:
        "{prefix}/partitions/covariance/{chromosome}/"
    run:
        "../scripts/calc_covariance_matrix.py"



