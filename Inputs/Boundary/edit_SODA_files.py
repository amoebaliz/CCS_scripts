import numpy as np
import netCDF4 as nc

# edit time
ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
dtot=0
newtime=np.zeros(len(ndays))

# establish sss_time values for cycling over 365d yr
# (i.e., middle of Jan, Feb, Mar ... Dec)
for nmon in range(len(ndays)):
    newtime[nmon] = dtot + ndays[nmon]/2.0
    dtot+=ndays[nmon]

ncfil = '/glade2/scratch2/edrenkar/CCS-inputs/CCS_bdry_soda3.4.1_1981-2010_clim.nc'
fid = nc.Dataset(ncfil,'a')

fid.variables['ocean_time'].units = 'Days'
fid.variables['ocean_time'].field = 'ocean_time, scalar, series'
fid.variables['ocean_time'][:] = newtime

fid.close()
