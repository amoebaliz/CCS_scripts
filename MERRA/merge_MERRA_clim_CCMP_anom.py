import numpy as np
import netCDF4 as nc
import datetime as dt


# loop over U and V files ?

dir = /glade2/scratch2/edrenkar/CCS-inputs/ 
clim_file = [drowned_MERRA_Uwind_1981-2010_MONTHLY_CLIM.nc, drowned_MERRA_Vwind_1981-2010_MONTHLY_CLIM.nc]
anom_file = anom_cat_1994.nc

fid_clim = nc.Dataset(clim_file)
fid_anom = nc.Dataset(anom_file)

## CREATE NEW NETCDF FILE ##
# DILEMA - what to read in for attributes








