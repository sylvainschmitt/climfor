rule get_tmf:
    input:
        "data/{area}/{area}.shp"
    output:
        "results/data/tmf_{area}.nc",
    log:
        "results/logs/get_tmf_{area}.log"
    benchmark:
        "results/benchmarks/get_tmf_{area}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    script:
      "../scripts/get_tmf.py"
