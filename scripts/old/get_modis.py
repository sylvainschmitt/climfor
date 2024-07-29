# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
file_ind = snakemake.output[0]
file_anom = snakemake.output[0]

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

# funs
def convert_to_df(ds):
    df = ds.to_dataframe()
    df.reset_index(inplace=True)
    return df
def calculate_year_month_cols(df):
    assert 'time' in df.columns, f"time should be in df.columns. Currently: {[c for c in df.columns]}"
    df['year'] = df.time.map(lambda x: x.year)
    df['month'] = df.time.map(lambda x: x.month)
    return df
def calculate_month_of_min_value(df, value_col):
    max_months = df.loc[df.groupby(["lat","lon","year"])[value_col].idxmin()]
    return max_months
def convert_dataframe_to_xarray(df, index_cols=['lat','lon']):
    out = df.set_index(index_cols).dropna()
    ds = out.to_xarray()
    return ds
def calculate_annual_month_of_min(ds, variable):
    df = convert_to_df(ds)
    df = calculate_year_month_cols(df)
    df = calculate_month_of_min_value(df, value_col=variable)
    ds_out = convert_dataframe_to_xarray(df, index_cols=['lat','lon','year'])
    return ds_out

# dryest month
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
dsp = xr.open_dataset("UCSB-CHG/CHIRPS/PENTAD", engine='ee')
dsp = dsp.sel(time=slice("2000-01-01", "2022-12-31"))
dsp = dsp.chunk(chunks = {'time':100, 'lat': 400, 'lon': 400})
dsp = dsp.transpose('time', 'lat', 'lon')
dsp = dsp.resample(time="1MS").sum()
regridder = xe.Regridder(dsp, ds_out, "bilinear")
dsp_r = regridder(dsp, keep_attrs=True)

# import matplotlib.pyplot as plt
# ds_test.sel(time="2001-01-01").tas.plot()
# plt.show()

# temperature
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
ds = xr.open_dataset("MODIS/061/MOD11A2", engine='ee')
ds = ds.sel(time=slice("2001-01-01", "2022-12-31"))
ds = ds.chunk(chunks = {'time':100, 'lat': 400, 'lon': 400})
ds = ds.transpose('time', 'lat', 'lon')
ds = ds.where(ds.QC_Day.isin([int("00000000",2), int("00000100",2)])) 
ds = ds.resample(time="1MS").mean(skipna = True)
regridder = xe.Regridder(ds, ds_out, "bilinear")
ds_r = regridder(ds, keep_attrs=True)

# join and filter
ds_all = xr.combine_by_coords([dsp_r, ds_r])
ds_all = ds_all.rename({"LST_Day_1km": "tas"})
ds_all = ds_all[["tas", "precipitation"]]
ds_max = calculate_annual_month_of_min(ds_all, variable='precipitation') 
ds_max = ds_max.sortby(['year', 'lat', 'lon'])
ds_max["tas"] = ds_max["tas"] * 0.02 - 273.15 # scale factor & degrees
ds_max[["tas"]].sel(year=2022).to_netcdf("results/data/modis_indices.nc")

# we excluded data where the estimated emissivity error was greater than 0.02 and where the LST error was greater than 1 K
    # Bits 0-1: Mandatory QA flags
    #     0: Pixel produced, good quality, not necessary to examine more detailed QA
    # Bits 2-3: Data quality flag
    #     0: Good data quality
    # Bits 4-5: Emissivity error flag
    #     0: Average emissivity error <= 0.01
    #     1: Average emissivity error <= 0.02
    # Bits 6-7: LST error flag
    #     0: Average LST error <= 1K
# we then selected the dryest month for each year and each pixel based on CHIRPS precipitation

# anomalies
ds_anom = ds_max.sel(year=slice(2020, 2022)).mean("year", skipna = True) -  ds_max.sel(year=slice(2001, 2003)).mean("year", skipna = True)
ds_anom[["tas"]].to_netcdf("results/data/modis_anomalies.nc")
