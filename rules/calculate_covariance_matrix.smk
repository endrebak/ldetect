rule effective_population_size:
    input:
        config["effective_population_sizes"]
    output:
        "{prefix}/partitions/{population}/samples.txt"
    run:
        df = pd.read_table(input[0], header=None)
        egs = df[df.Population.isin(pop_as_list())].EGS.mean()

        with open(output[0], "w+") as egs_h:
            egs_h.write(str(egs) + "\n")


checkpoint tabix_intervals:
    input:
        genetic_maps = rules.interpolate_genetic_maps.output,
        variants = rules.fetch_variants.output
    output:
        "{prefix}/partitions/{population}/{chromosome}/{intervals}.vcf.gz"
    run:
        df = pd.read_table(input.genetic_maps)
        w = wildcards

        prefix = w.prefix
        population = w.population
        chromosome = w.chromosome

        c = chromosome.replace("chr", "")
        for _, (start, end) in df.iterrows():
            o = f"{prefix}/partitions/{population}/{chromosome}/{start}_{end}.vcf.gz"
            cmd = f"tabix {input.variants} {c}:{start}-{end} | gzip > {o}"
            print(cmd)
            shell(cmd)



rule calculate_covariance_matrix:
    input:
        genetic_map = rules.interpolate_genetic_maps,
        samples = rules.number_individuals_in_reference_panel.output.samples,
        effective_population_size = rules.effective_population_size.output
        # interval = checkpoint.
    output:
        "{prefix}/partitions/covariance/samples.txt"
    script:
        "../scripts/calc_covariance_matrix.py"



rule matrix_to_vector:
    input:
        glob_intervals
    output:
        


def glob_intervals(w):

    checkpoint_output = checkpoints.tabix_intervals.get(**w).output[0]
    file_template = "{prefix}/partitions/{population}/{chromosome}/{{intervals}}.vcf.gz".format(**w)
    intervals = glob_wildcards(file_template)

    return expand(file_template, intervals=intervals)


