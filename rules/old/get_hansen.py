rule get_hansen:
    input:
        "data/{area}/{area}.shp"
    output:
        "results/data/hansen_{area}.nc",
    log:
        "results/logs/get_hansen_{area}.log"
    benchmark:
        "results/benchmarks/get_hansen_{area}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    script:
      "../scripts/get_hansen.py"
