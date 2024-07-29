# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
amazonia = snakemake.input
nc = snakemake.output

# test
# amazonia = "results/data/amazonia"
# amazonia = "results/data/guaviare"
# amazonia = "results/data/capricho"
# nc = "results/data/tmf/tmf.nc"
        
# libs
import geopandas as gp
import ee
import xarray as xr
import numpy as np

# code
area = gp.read_file(amazonia).dissolve()
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
leg = ee.Geometry.Rectangle(area.bounds.minx[0], area.bounds.miny[0], area.bounds.maxx[0], area.bounds.maxy[0])
ic = ee.ImageCollection("projects/JRC/TMF/v1_2022/AnnualChanges")
ds = xr.open_mfdataset(
        [ic],
        engine='ee',
        projection=ic.first().select(0).projection(),
        geometry=leg
).sel(time=2).drop('time')
ds2 = ds.to_array("year", name="tmf").to_dataset().assign_coords(
    year=np.arange(1990,2023)
).chunk({"lon": 1000, "lat": 1000, "year": 10})
# ds2 = ds2.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
# ds2 = ds2.rio.write_crs("epsg:4362")
# ds2 = ds2.rio.clip(area.geometry.values, area.crs)
ds2.to_netcdf(nc)
# val 1. Undisturbed Tropical moist forest (TMF)  
# val 2. Degraded TMF  
# val 3. Deforested land  
# val 4. Forest regrowth  
# val 5. Permanent or seasonal Water  
# val 6. Other land cover  

# plot
# from matplotlib import pyplot as plt
# area.plot()
# ds2.sel(year=1990)["tmf"].plot()
# plt.show()
