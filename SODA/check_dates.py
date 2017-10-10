# Auxiliary sanity check of the relative temporal position
# of the monthly climatologies. i.e., does the jan 81-10 clim 
# fall before the feb 81-10 clim; the actual dates are meaningless


import numpy as np
import netCDF4 as nc
import datetime as dt

ref = dt.datetime(1980,01,01)

#for yr in range(1980,2009+1):
#    print yr
#    for mon in range(1,12+1):
#        filstr = './' + str(yr) + '/soda3.3.1_mon_' + str(yr) + '_' + str(mon).zfill(2) + '.nc'
#        fid = nc.Dataset(filstr)

#        time = fid.variables['time'][:]
#        print ref + dt.timedelta(days=time[0])

for mon in range(1,12+1):
    print mon
    time = np.zeros(31)
    for yr in range(1980,2010+1):
        filstr = '../' + str(yr) + '/soda3.3.1_mon_' + str(yr) + '_' + str(mon).zfill(2) + '.nc'
        fid = nc.Dataset(filstr)

        time[yr-1980] = fid.variables['time'][:]
        print ref + dt.timedelta(days=fid.variables['time'][0])
    print 'MEAN', ref + dt.timedelta(days=np.mean(time))






