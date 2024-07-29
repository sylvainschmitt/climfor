rule get_modis:
    input:
        "results/limits/limits.shp"
    output:
        "results/data/modis.nc"
    log:
        "results/logs/get_modis.log"
    benchmark:
        "results/benchmarks/get_modis.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    script:
      "../scripts/get_modis.py"
