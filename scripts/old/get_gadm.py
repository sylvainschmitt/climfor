# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
area = snakemake.params.area
area_file = snakemake.output[0]

# test
# area="Guainía"

# libs
import pygadm

# code
code = pygadm.AdmNames(area).GID_1[0]
gdf = pygadm.AdmItems(admin = code)
gdf.to_file(area_file)
