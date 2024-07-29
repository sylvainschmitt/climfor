# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
amazonia = snakemake.input
nc = snakemake.output

# test
# amazonia = "results/data/amazonia"
# amazonia = "results/data/guaviare"
amazonia = "results/data/capricho"
nc = "results/data/precipitation/chirps.nc"
        
# libs
import geopandas as gp
import ee
import xarray as xr
        
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
ds.pr.attrs = {'standard_name': 'precipitation', 
                'long_name': 'Monthly precipitation',
                'units': 'mm month-1', 
                'explanation' : 'Precipitation in the earth\'s atmosphere, monthly means precipitation of water in all phases.'}
ds.to_netcdf(nc)

# plot
# from matplotlib import pyplot as plt
# area.plot()
# ds.sel(time="2001-01-01")["pr"].plot()
# plt.show()

