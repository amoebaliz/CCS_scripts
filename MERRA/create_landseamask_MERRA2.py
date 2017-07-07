import os
import numpy as np
import netCDF4 as nc

ncfil = 'MERRA2_101.const_2d_asm_Nx.00000000.nc4'
os_commd = 'wget goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2C0NXASM.5.12.4/1980/' + ncfil 
#os.system(os_commd)

f = nc.Dataset(ncfil)

frac_ocean   = f.variables['FROCEAN'][:].squeeze()

lon  = f.variables['lon'][:]
lat  = f.variables['lat'][:]

f.close()

# ------------------------------ Create binary land sea mask --------------------------------------------------
ny, nx = frac_ocean.shape
lsm = np.zeros((ny,nx))

threshold = 0.995
lsm[np.where(frac_ocean >= threshold)] = 1

# ------------------------------ Write mask to netcdf ---------------------------------------------------------
fout = nc.Dataset('./MERRA2_lsm.nc','w',format='NETCDF3_CLASSIC')
fout.description = 'MERRA land sea mask '

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
