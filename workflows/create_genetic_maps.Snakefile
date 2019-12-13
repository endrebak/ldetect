prefix = config["prefix"]

sample_sheet_path = config["sample_sheet"]

rule download_maps:
    output:
        "{prefix}/genetic_map/{population}.gz"
    shell:
        "curl {}"
