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
nc = "results/data/precipitation/chirps.nc"
        
# libs
import geopandas as gp
import ee
import xarray as xr
import numpy as np
from scipy import stats 
        
# daily pr
area = gp.read_file(amazonia).dissolve()
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
ic = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
leg1 = ee.Geometry.Rectangle(area.bounds.minx[0], area.bounds.miny[0], area.bounds.maxx[0], area.bounds.maxy[0])
ds = xr.open_mfdataset(
        [ic],
        engine='ee',
        projection=ic.first().select(0).projection(),
        geometry=leg1
).chunk({"lon": 1000, "lat": 1000, "time": 10})
ds = ds.transpose('time', 'lat', 'lon')
ds = ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
ds = ds.rename({'precipitation' : 'pr'})
ds.to_netcdf(nc)

# yearly pr
ds_year = ds.resample(time="Y").sum()
ds_year["time"] = ds_year.time.dt.year
ds_year = ds_year.rename({"time": "year"})
ds_year.to_netcdf("results/data/precipitation/chirps_year.nc")

# yearly anom
ds_year_anom = ds_year
ds_year_anom['pr'] = ((ds_year.pr - ds_year.mean(dim="year").pr)+0.1)/(ds_year.std(dim="year").pr+0.1)
ds_year_anom = ds_year_anom.rename({'pr' : 'pr_stdanom'})
ds_year_anom.to_netcdf("results/data/precipitation/chirps_anom.nc")

# yearly tau
ds_year_anom = ds_year_anom.chunk({"lon": 1000, "lat": 1000, "year": -1}) # to work across years
# ds_year_anom = ds_year_anom.compute()

def tau(x, y):
        tau, p_value = stats.kendalltau(x, y)
        return np.array([tau, p_value])

stats = xr.apply_ufunc(tau, 
                       ds_year_anom['pr_stdanom'], 
                       ds_year_anom['year'], 
                       input_core_dims=[['year'], ['year']],
                       output_core_dims=[["parameter"]],
                       vectorize=True,
                       dask="parallelized",
                       output_dtypes=['float64'],
                       output_sizes={"parameter": 2},
                      )

ds_year_tau = ds_year_anom[['lon', 'lat']]
ds_year_tau["tau"] = stats[:,:,0]
ds_year_tau["pval"] = stats[:,:,1]
ds_year_tau.to_netcdf("results/data/precipitation/chirps_tau.nc")

# plot
from matplotlib import pyplot as plt
ds_year_tau["pval"].plot()
plt.show()
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
nc = "results/data/precipitation/chirps.nc"
        
# libs
import geopandas as gp
import ee
import xarray as xr
import numpy as np
from scipy import stats 
        
# daily pr
area = gp.read_file(amazonia).dissolve()
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
ic = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
leg1 = ee.Geometry.Rectangle(area.bounds.minx[0], area.bounds.miny[0], area.bounds.maxx[0], area.bounds.maxy[0])
ds = xr.open_mfdataset(
        [ic],
        engine='ee',
        projection=ic.first().select(0).projection(),
        geometry=leg1
).chunk({"lon": 1000, "lat": 1000, "time": 10})
ds = ds.transpose('time', 'lat', 'lon')
ds = ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
ds = ds.rename({'precipitation' : 'pr'})
ds.to_netcdf(nc)

# yearly pr
ds_year = ds.resample(time="Y").sum()
ds_year["time"] = ds_year.time.dt.year
ds_year = ds_year.rename({"time": "year"})
ds_year.to_netcdf("results/data/precipitation/chirps_year.nc")

# yearly anom
ds_year_anom = ds_year
ds_year_anom['pr'] = ((ds_year.pr - ds_year.mean(dim="year").pr)+0.1)/(ds_year.std(dim="year").pr+0.1)
ds_year_anom = ds_year_anom.rename({'pr' : 'pr_stdanom'})
ds_year_anom.to_netcdf("results/data/precipitation/chirps_anom.nc")

# yearly tau
ds_year_anom = ds_year_anom.chunk({"lon": 1000, "lat": 1000, "year": -1}) # to work across years
# ds_year_anom = ds_year_anom.compute()

def tau(x, y):
        tau, p_value = stats.kendalltau(x, y)
        return np.array([tau, p_value])

stats = xr.apply_ufunc(tau, 
                       ds_year_anom['pr_stdanom'], 
                       ds_year_anom['year'], 
                       input_core_dims=[['year'], ['year']],
                       output_core_dims=[["parameter"]],
                       vectorize=True,
                       dask="parallelized",
                       output_dtypes=['float64'],
                       output_sizes={"parameter": 2},
                      )

ds_year_tau = ds_year_anom[['lon', 'lat']]
ds_year_tau["tau"] = stats[:,:,0]
ds_year_tau["pval"] = stats[:,:,1]
ds_year_tau.to_netcdf("results/data/precipitation/chirps_tau.nc")

# plot
from matplotlib import pyplot as plt
ds_year_tau["pval"].plot()
plt.show()
