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
        ds = ds.where(ds.b1 > 0)
        ds['b1'] = ds.b1*0.02 - 273.15  # K to °C
        if(var == "tas"):
                std_name = "temperature"
        if(var == "tasmin"):
                std_name = "minimum temperature"
        if(var == "tasmax"):
                std_name = "maximum temperature"
        if(var == "tas"):
                type = "mean"
        if(var == "tasmin"):
                type = "minimum"
        if(var == "tasmax"):
                type = "maximum"
        ds.b1.attrs = {'standard_name': std_name + ' at surface', 
                        'long_name': 'Monthly ' + type + ' daily air temperature',
                        'units': '°C', 
                        'explanation' : 'Monthly ' + type + ' air temperatures at 2 meters.'}
        ds = ds.rename({'b1': var})
        return(ds)

def get_gshtd(bounds):
        leg1 = ee.Geometry.Rectangle(bounds.minx[0], bounds.miny[0], bounds.maxx[0], bounds.maxy[0])
        ds_tas = get_var("tas", leg1)
        ds_tasmin = get_var("tasmin", leg1)
        ds_tasmax = get_var("tasmax", leg1)
        ds = xr.merge([ds_tas, ds_tasmax, ds_tasmin])
        return(ds)

# code
area = gp.read_file(amazonia).dissolve()
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
gshtd = get_gshtd(area.bounds).chunk({"lon": 1000, "lat": 1000, "time": 10})
gshtd.to_zarr(zarr)

# check
# ds = xr.open_dataset(zarr)
# from matplotlib import pyplot as plt
# ds.sel(time="2001-01-01").pr.plot()
# plt.show()
