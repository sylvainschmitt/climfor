# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
amazonia = snakemake.input
nc = snakemake.output

# test
amazonia = "results/data/amazonia"
zarr = "results/data/temperature/gshtd.zarr"
        
# libs
import geopandas as gp
import ee
import xarray as xr

# funs
def get_var(var, leg):
        if(var == "tas"):
                col = "TMEAN"
        if(var == "tasmin"):
                col = "TMIN"
        if(var == "tasmax"):
                col = "TMAX"
        ic = ee.ImageCollection("projects/sat-io/open-datasets/GSHTD/" + col)
        ds = xr.open_dataset(
                ic,
                engine='ee',
                projection=ic.first().select(0).projection(),
                geometry=leg
        )
        ds = ds.transpose('time', 'lat', 'lon')
        ds = ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
        ds = ds.rio.write_crs("epsg:4362")
        return(ds)

# each raw
area = gp.read_file(amazonia).dissolve()
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com', project="silver-link-350720")
leg1 = ee.Geometry.Rectangle(area.bounds.minx[0], area.bounds.miny[0], area.bounds.maxx[0], area.bounds.maxy[0])
ds_tas = get_var("tas", leg1).chunk({"lon": 1000, "lat": 1000, "time": 10})
ds_tas.to_zarr("results/data/temperature/gshtd_tas.zarr")
ds_tasmin = get_var("tasmin", leg1).chunk({"lon": 1000, "lat": 1000, "time": 10})
ds_tasmin.to_zarr("results/data/temperature/gshtd_tasmin.zarr")
ds_tasmax = get_var("tasmax", leg1).chunk({"lon": 1000, "lat": 1000, "time": 10})
ds_tasmax.to_zarr("results/data/temperature/gshtd_tasmax.zarr")

# assemble & prep each (order to be defined), save as a single zarr and remove the tmp
# dask_jobqueue could be used to increase speed
import xarray as xr
ds_tas = xr.open_zarr("results/data/temperature/gshtd_tas.zarr")
ds_tas.b1.attrs = {'standard_name': "temperature" + ' at surface', 
                'long_name': 'Monthly ' + "mean" + ' daily air temperature',
                'units': '°C', 
                'explanation' : 'Monthly ' + "mean" + ' air temperatures at 2 meters.'}
ds_tas = ds_tas.rename({'b1': "tas"})

ds_tasmin = xr.open_zarr("results/data/temperature/gshtd_tasmin.zarr")
ds_tasmin.b1.attrs = {'standard_name': " minimum temperature" + ' at surface', 
                'long_name': 'Monthly ' + "minimum" + ' daily air temperature',
                'units': '°C', 
                'explanation' : 'Monthly ' + "minimum" + ' air temperatures at 2 meters.'}
ds_tasmin = ds_tasmin.rename({'b1': "tasmin"})

ds_tasmax = xr.open_zarr("results/data/temperature/gshtd_tasmax.zarr")
ds_tasmax.b1.attrs = {'standard_name': "maximum temperature" + ' at surface', 
                'long_name': 'Monthly ' + "maximum" + ' daily air temperature',
                'units': '°C', 
                'explanation' : 'Monthly ' + "maximum" + ' air temperatures at 2 meters.'}
ds_tasmax = ds_tasmax.rename({'b1': "tasmax"})

ds = xr.merge([ds_tas, ds_tasmax, ds_tasmin])
ds = ds.chunk({"lon": 1000, "lat": 1000, "time": 10})
ds.to_zarr("results/data/temperature/gshtd_raw.zarr")

# remove null and rescale raw, save as a single zarr and remove the tmp
# dask_jobqueue could be used to increase speed
import xarray as xr
ds = xr.open_zarr("results/data/temperature/gshtd_raw.zarr")
ds = ds.where(ds.tas > 0)
ds['tas'] = ds.tas*0.02 - 273.15  # K to °C
ds['tasmax'] = ds.tas*0.02 - 273.15  # K to °C
ds['tasmin'] = ds.tas*0.02 - 273.15  # K to °C
ds.to_zarr("results/data/temperature/gshtd.zarr")

# check
# import xarray as xr
# from matplotlib import pyplot as plt
# ds = xr.open_zarr("results/data/temperature/gshtd.zarr")
# ds.sel(time="2001-01-01").tas.plot()
# plt.savefig("test.png")
# plt.show()
