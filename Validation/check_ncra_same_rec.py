import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

ncfil1 = 'SST_am_01_clim.nc' 
ncfil2 = 'SST_pm_01_clim.nc' 
ncfil3 = 'SST_all_01_clim.nc'

sst1 = nc.Dataset(ncfil1).variables['BSST'][:].squeeze()
sst2 = nc.Dataset(ncfil2).variables['BSST'][:].squeeze()
sst3 = nc.Dataset(ncfil3).variables['BSST'][:].squeeze()

print (sst1+sst2)/2. == sst3


