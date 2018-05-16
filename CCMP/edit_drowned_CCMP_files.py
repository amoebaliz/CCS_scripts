import netCDF4 as nc

filvar = ['u10', 'v10']
invar = ['u_anom','v_anom']
invar2 = ['Uwind', 'Vwind']

for nvar in range(len(invar)):

    ncfil1 = '/Users/elizabethdrenkard/external_data/ERAinterim/drowned/drowned_ERAi_' + filvar[nvar] + '_1981-2010_monthly_clim.nc'
    ncfil2 = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/CCMP/drowned_CCMP_' + filvar[nvar][0] + '_anom.nc'

    fid1 = nc.Dataset(ncfil1)
    fid2 = nc.Dataset(ncfil2,'a')

    fid2.variables['time'].units = fid1.variables['time'].units    
    fid2.variables['time'].cycle_length = fid1.variables['time'].cycle_length 

    fid2.variables['lon'].long_name = fid1.variables['lon'].long_name
    fid2.variables['lat'].long_name = fid1.variables['lat'].long_name

    fid2.variables[invar[nvar]].missing_value = fid1.variables[invar2[nvar]].missing_value
    fid2.variables[invar[nvar]].time = 'time'

    fid1.close()
    fid2.close()
