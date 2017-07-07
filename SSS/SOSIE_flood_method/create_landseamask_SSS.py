import numpy as np
import netCDF4 as nc

ncfil = 'phc3.0_monthly.nc'
f = nc.Dataset(ncfil)

salt   = f.variables['salt'][0,:2,:,:].squeeze()
salt   = np.mean(salt,axis=0)

lon  = f.variables['lon'][:]
lat  = f.variables['lat'][:]

f.close()
# ------------------------------ Create binary land sea mask --------------------------------------------------
ny, nx = salt.shape
lsm = np.zeros((ny,nx))

threshold = 50
lsm[np.where( salt <= threshold)] = 1
# ------------------------------ Write mask to netcdf ---------------------------------------------------------
fout = nc.Dataset('./phc3.0_lsm.nc','w',format='NETCDF3_CLASSIC')
fout.description = 'PHC land sea mask '

# dimensions
fout.createDimension('lat', ny)
fout.createDimension('lon', nx)

# variables
latitudes  = fout.createVariable('lat', 'f8', ('lat',))
longitudes = fout.createVariable('lon', 'f8', ('lon',))
mask       = fout.createVariable('lsm', 'i4', ('lat','lon',))

# data
latitudes[:]    = lat
longitudes[:]   = lon
mask[:,:]       = lsm

# close
fout.close()
