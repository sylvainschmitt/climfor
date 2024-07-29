rule crop_chirps:
    input:
        "data/{area}/{area}.shp",
        "results/data/chirps.nc"
    output:
        "results/data/chirps_{area}.nc",
    log:
        "results/logs/crop_chirps_{area}.log"
    benchmark:
        "results/benchmarks/crop_chirps_{area}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    script:
      "../scripts/crop_chirps.py"
