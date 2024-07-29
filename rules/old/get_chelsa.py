rule get_chelsa:
    input:
        "data/{area}/{area}.shp"
    output:
        "results/data/chelsa_{area}.nc",
    log:
        "results/logs/get_chelsa_{area}.log"
    benchmark:
        "results/benchmarks/get_chelsa_{area}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    script:
      "../scripts/get_chelsa.py"
