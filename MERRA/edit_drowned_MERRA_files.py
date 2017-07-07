import netCDF4 as nc

#invar = ['Qair','Pair','Tair','Uwind','Vwind','swrad','lwrad_down','rain']
invar = ['Uwind','Vwind']
for nvar in range(len(invar)):
   # ncfil1 = 'MERRA_' + invar[nvar] + '_1981-2010_MONTHLY_CLIM.nc'
   # ncfil2 = './drowned/drowned_MERRA_' + invar[nvar] + '_1981-2010_MONTHLY_CLIM.nc'
    ncfil1 = 'MERRA_' + invar[nvar] + '_10m_1981-2010_MONTHLY_CLIM.nc'
    ncfil2 = './drowned/drowned_MERRA_' + invar[nvar] + '_10M_1981-2010_MONTHLY_CLIM.nc'
    fid1 = nc.Dataset(ncfil1)
    fid2 = nc.Dataset(ncfil2,'a')

    fid2.variables['time'].units = fid1.variables['time'].units    
    fid2.variables['time'].cycle_length = fid1.variables['time'].cycle_length 

    fid2.variables['lon'].long_name = fid1.variables['lon'].long_name
    fid2.variables['lat'].long_name = fid1.variables['lat'].long_name

    fid2.variables[invar[nvar]].missing_value = fid1.variables[invar[nvar]].missing_value
    fid2.variables[invar[nvar]].time = 'time'

    fid1.close()
    fid2.close()
