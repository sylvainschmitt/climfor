# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test

# libs
import ee
import xarray as xr
import xesmf as xe
import numpy as np

# grid
ds_out = xr.Dataset(
    {
        "lat": (["lat"], np.arange(-10, 10, 0.05), {"units": "degrees_north"}),
        "lon": (["lon"], np.arange(-80, -60, 0.05), {"units": "degrees_east"}),
    }
)

# get
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
ds = xr.open_dataset("UCSB-CHG/CHIRPS/PENTAD", engine='ee')
ds = ds.sel(time=slice("2000-01-01", "2022-12-31"))
ds = ds.chunk(chunks = {'time':100, 'lat': 400, 'lon': 400})
ds = ds.transpose('time', 'lat', 'lon')
# import matplotlib.pyplot as plt
# ds.sel(time="2001-01-01").precipitation.plot()
# plt.show()

# indices
ds = ds.groupby("time.year").sum()
regridder = xe.Regridder(ds, ds_out, "bilinear")
ds_r = regridder(ds, keep_attrs=True)
ds_r.to_netcdf("results/data/chirps_indices.nc")

# anomalies
ds_anom = ds.sel(year=slice(2020, 2022)).mean("year") -  ds.sel(year=slice(2000, 2002)).mean("year")
regridder = xe.Regridder(ds_anom, ds_out, "bilinear")
ds_r = regridder(ds_anom, keep_attrs=True)
ds_r.to_netcdf("results/data/chirps_anomalies.nc")
