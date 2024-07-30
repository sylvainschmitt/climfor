# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
amazonia = snakemake.input
nc = snakemake.output

# test
# amazonia = "results/data/amazonia"
amazonia = "results/data/guaviare"
# amazonia = "results/data/capricho"
# nc = "results/data/tmf/tmf.nc"
nc = "results/data/forest/deforestation_year.nc"
tif = "results/data/forest/deforestation_year.tif"


# libs
import geopandas as gp
import ee
import xarray as xr
import numpy as np
import xesmf as xe

# code
area = gp.read_file(amazonia).dissolve()
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
leg = ee.Geometry.Rectangle(area.bounds.minx[0], area.bounds.miny[0], area.bounds.maxx[0], area.bounds.maxy[0])
ic = ee.ImageCollection("projects/JRC/TMF/v1_2022/DeforestationYear")
ds = xr.open_mfdataset(
        [ic],
        engine='ee',
        projection=ic.first().select(0).projection(),
        geometry=leg
).sel(time=2).drop('time').rename({"constant" : "year"})
ds = ds.transpose('lat', 'lon')
ds = ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
ds = ds.rio.write_crs("epsg:4362")
ds = ds.where(ds.year > 1989)
ds.rio.to_raster(tif)   
