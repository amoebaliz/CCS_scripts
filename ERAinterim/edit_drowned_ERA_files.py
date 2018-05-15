import netCDF4 as nc

#filvar = ['q2','msl','t2','u10','v10','radsw','radlw','precip']
#invar = ['Qair','Pair','Tair','Uwind','Vwind','swrad','lwrad_down','rain']
filvar = ['radsw']
invar = ['swrad']

#invar = ['Uwind','Vwind']
for nvar in range(len(invar)):
    #ncfil1 = '/Users/elizabethdrenkard/external_data/ERAinterim/global/' + filvar[nvar] + '_1981-2010_monthly_clim.nc'
    #ncfil2 = '/Users/elizabethdrenkard/external_data/ERAinterim/global/drowned/drowned_ERAi_' + filvar[nvar] + '_1981-2010_monthly_clim.nc'
    #ncfil1 = '/glade/p/work/edrenkar/external_data/ERAinterim/climatologies/global/'+ filvar[nvar] + '_1981-2010_monthly_clim.nc'
    #ncfil2 = '/glade/p/work/edrenkar/external_data/ERAinterim/climatologies/global/drowned/drowned_ERAi_' + filvar[nvar] + '_1981-2010_monthly_clim.nc'
    ncfil1 = '/Users/liz.drenkard/external_data/ERAinterim/Forecast/climatologies/'+ filvar[nvar] + '_1981-2010_monthly_clim.nc'
    ncfil2 = '/Users/liz.drenkard/external_data/ERAinterim/Forecast/climatologies/drowned/drowned_ERAi_' + filvar[nvar] + '_1981-2010_monthly_clim.nc'
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
