rule get_amazonia:
    output:
      directory("results/data/amazonia/"),
      temp("results/data/amazonia/amazonia.zip")
    log:
      "results/logs/get_amazonia.log"
    benchmark:
      "results/benchmarks/get_amazonia.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    params:
      amazonia=config["amazonia"]
    shell:
      "wget {params.amazonia} -O {output[1]} ; "
      "unzip {output[1]} -d {output[0]}"
