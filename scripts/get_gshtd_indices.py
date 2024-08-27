# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
amazonia = snakemake.input
nc = snakemake.output

# test
zarr = "results/data/temperature/gshtd.zarr"
        
# libs
import xarray as xr
import numpy as np
from scipy import stats 
        
# daily pr
ds = xr.open_zarr(zarr)

# yearly pr
ds_year = ds.resample(time="Y").mean()
ds_year["time"] = ds_year.time.dt.year
ds_year = ds_year.rename({"time": "year"})
ds_year.to_netcdf("results/data/temperature/gshtd_year.nc")

# yearly anom
ds_year_anom = ds_year
ds_year_anom = ds_year - ds_year.mean(dim="year")
ds_year_anom.to_netcdf("results/data/temperature/gshtd_anom.nc")

# yearly stdanom
ds_year_anom = ds_year
ds_year_anom = (ds_year - ds_year.mean(dim="year"))/ds_year.std(dim="year")
ds_year_anom.to_netcdf("results/data/temperature/gshtd_stdanom.nc")

# yearly tau
ds_year_anom = ds_year_anom.chunk({"lon": 1000, "lat": 1000, "year": -1}) # to work across years
# ds_year_anom = ds_year_anom.compute()

import xarray as xr
import numpy as np
from scipy import stats 

# ds_year_anom = xr.open_dataset("results/data/precipitation/.nc")

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
# ds_year_tau.to_netcdf("results/data/precipitation/.nc")
