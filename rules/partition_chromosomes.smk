rule number_individuals_in_reference_panel:
    input:
        config["sample_info"]
    output:
        n_ind = "{prefix}/partitions/{population}/n_individuals.txt",
        samples = "{prefix}/partitions/{population}/samples.txt"
    run:
        df = pd.read_table(input[0])

        pop = df.Population.isin(pop_as_list())

        pop = len(pop.drop_duplicates("Sample"))

        with open(output["n_ind"], "w+") as o:
            o.write(str(pop) + "\n")

        pop.Sample.to_frame().to_csv(output["samples"], index=False, sep="\t", header=False)


rule partition_chromosomes:
    input:
        genetic_map = rules.interpolate_genetic_maps.output,
        number_individuals = rules.number_individuals_in_reference_panel.output.n_ind
    output:
        "{prefix}/partitions/{population}/{chromosome}.gz"
    script:
        "scripts/partition_chromosome.py"


