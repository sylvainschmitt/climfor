rule get_gadm:
    output:
        directory("results/area/{area}")
    log:
        "results/logs/get_gadm_{area}.log"
    benchmark:
        "results/benchmarks/get_gadm_{area}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/gadm.yml"
    params:
        area="{area}"
    script:
      "../scripts/get_gadm.py"
