rule get_chirps:
    input:
        "results/limits/limits.shp"
    output:
        "results/data/chirps.nc"
    log:
        "results/logs/get_chirps.log"
    benchmark:
        "results/benchmarks/get_chirps.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    script:
      "../scripts/get_chirps.py"

