configfile: "config/config.yml"

# rules #
rule all:
   input:
      # Amazonia
      "results/data/amazonia/"

## rules ##
include: "rules/get_amazonia.py"
