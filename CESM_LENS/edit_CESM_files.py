import matplotlib
matplotlib.use('Agg')
import numpy as np
import netCDF4 as nc
from datetime import datetime
import sys

#filvar = ['Q','PSL','TREFHT','U','V','FSNS','FLDS','PREC_F']
#filvar = ['q2','msl','t2','u10','v10','radsw','radlw','precip']

ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
dtot=0
dtim = np.zeros(12)
for nt in range(12):
    dtim[nt] = dtot + ndays[nt]/2.0
    dtot+=ndays[nt]

for nvar in range(len(filvar)):
           
        # edit ROMS forcing file
        #ncfile = 'ERAi_' + filvar[nvar] + '_1981-2010_monthly_clim.nc'
        
        fid    = nc.Dataset(ncfile, 'a', format='NETCDF3_CLASSIC')
        fid.variables['time'].units = 'Days'
        fid.variables['time'].calendar = 'NOLEAP'
        fid.variables['time'][:] = dtim
        fid.close()
