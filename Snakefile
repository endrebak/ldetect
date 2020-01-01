import pandas as pd


def pop_as_list():
    pop = config["population"]
    if isinstance(pop, str):
        return [pop]
    elif isinstance(pop, list):
        return pop
    else:
        values = []
        for v in list(pop.values())[0]:
            values.append(v)

        return values


sample_sheet_path = config["sample_sheet"]
prefix = config["prefix"]

ss = pd.read_table(sample_sheet_path)

si = pd.read_table("data/sample_info.tsv")


def get_number_of_individuals_in_reference_panel(w):

    return len(si[si.Population == w.population])


chromosomes = ["chr" + str(i) for i in range(1, 23)]
chromosomes = ["chr1"] #, "chr22"]
prefix = config["prefix"]
populations = pop_as_list()
# populations = pd.read_table(config["recombination_rates"], header=None, squeeze=True).to_list()

regex = lambda l: "|".join([str(w) for w in l])


wildcard_constraints:
    chromosome = regex(chromosomes),
    prefix = prefix,
    population = regex(pop_as_list())


for rule in ["interpolate_genetic_maps", "partition_chromosomes", "calculate_covariance_matrix", "find_minima"]:
    # print(rule)
    include: "rules/" + rule + ".smk"


def aexpand(f):

    if not isinstance(f, str):
        f = str(f)

    import re
    get_wildcards = r"\{(.*?)\}"
    ws = list(res.group(1) for res in re.finditer(get_wildcards, f))

    _globals = globals()
    keys_found = {w: w for w in ws if w in _globals}
    keys_found_with_s = {w: w + "s" for w in ws if w + "s" in _globals}

    keys_found_with_s.update(keys_found)
    keys_found = keys_found_with_s

    k_v = {k: _globals[keys_found[k]] for k in keys_found}

    return expand(f, **k_v)

rule all:
    input:
        # aexpand(rules.partition_chromosomes.output)
        # aexpand(checkpoints.covariance_matrix.output)
        # aexpand("{prefix}/partition_covariances/{population}/{chromosome}/")
        # aexpand("{prefix}/partitions/covariance/{population}/{chromosome}.gz")
        aexpand("{prefix}/partitions/covariance_much_ram/{population}/{chromosome}.txt")
        # aexpand("{prefix}/minima/{population}/{chromosome}.gz")
        # aexpand("{prefix}/collected_covariances/{population}/{chromosome}.gz")
        # aexpand(rules.fetch_variants.output)

# rule download_chromosomes_vcf:
#     output:
#         "{prefix}/download/{chromosome}.gz"
#     params:
#         url = lambda w: ss[ss.Chromosome == w.chromosome].URL.iloc[0]
#     shell:
#         "curl {params.url} > {output[0]}"


# rule liftover_recombination_maps:
#     input:
#         config["genetic_map"]
#     output:
#         "{prefix}/download/{chromosome}.gz"



