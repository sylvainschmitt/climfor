rule get_area:
    input:
        expand("results/area/{area}", area=config["regions"])
    output:
        directory("results/limits"),
        "results/limits/limits.shp"
    log:
        "results/logs/get_area.log"
    benchmark:
        "results/benchmarks/get_area.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    script:
      "../scripts/get_area.py"
