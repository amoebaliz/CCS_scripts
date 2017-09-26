# open each file
# change the time value
# change the time unit

import numpy as np
import netCDF4 as nc

ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
time_vals = np.zeros(12)
dtot = 0

for n in range(12):
    time_vals[n] = ndays[n]/2.0 + dtot
    dtot += ndays[n]

dir_an = '/Users/liz.drenkard/external_data/ERAinterim/Forecast/climatologies/'
era_vars = ['t2','q2','u10','v10','msl','radsw','radlw','precip']

for var in era_vars:
    ncfil = dir_an + var + '_1981-2010_monthly_clim.nc'
    fid = nc.Dataset(ncfil,'a')
    fid.variables['time'].units = 'Days'
    fid.variables['time'].calendar = 'NOLEAP'
    fid.variables['time'][:] = time_vals    
    fid.close()
