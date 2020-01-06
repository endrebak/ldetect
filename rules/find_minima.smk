

rule find_minima:
    input:
        intervals = get_intervals,
        vector = rules.matrix_to_vector_much_ram.output[0]
    output:
        "{prefix}/minima/{population}/{chromosome}.gz"
    shell:
        "python scripts/find_minima.py {input.vector} {input.intervals}"


