# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
areas = snakemake.input
file = snakemake.output[0]

# code
import geopandas as gp
import pandas as pd
pd.concat(list(map(gp.read_file, areas))).to_file(file)
