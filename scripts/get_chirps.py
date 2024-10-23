# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
amazonia = snakemake.input
nc = snakemake.output

# test
amazonia = "results/data/amazonia"
zarr = "results/data/precipitation/chirps.zarr"
        
# libs
import geopandas as gp
import ee
import xarray as xr
        
# daily pr
area = gp.read_file(amazonia).dissolve()
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
ic = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
leg1 = ee.Geometry.Rectangle(area.bounds.minx[0], area.bounds.miny[0], area.bounds.maxx[0], area.bounds.maxy[0])
ds = xr.open_dataset(
        ic,
        engine='ee',
        projection=ic.first().select(0).projection(),
        geometry=leg1
).chunk({"lon": 1000, "lat": 1000, "time": 10})
ds = ds.transpose('time', 'lat', 'lon')
ds = ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
ds = ds.rename({'precipitation' : 'pr'})
ds.to_zarr(zarr)

# check
# ds = xr.open_dataset(zarr)
# from matplotlib import pyplot as plt
# ds.sel(time="1981-01-01").pr.plot()
# plt.show()

ds = xr.open_dataset(zarr)
import numpy as np
x_coord = -54.9588900
y_coord = -2.85667
ds2 = ds.interp({'lon':x_coord, 'lat':y_coord}, method='linear')
ds2.to_netcdf("tapajos_chirps.nc")

val = np.diagonal(ds.interp({'lon':x_coord, 'lat':y_coord}, method='linear').values[0,:,:])
print(val)
