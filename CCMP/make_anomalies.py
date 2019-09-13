import numpy as np
import netCDF4 as nc
import datetime as dt

# LOAD CLIMATOLOGY FILE
#fidclim = nc.Dataset('CCMP_1988-2010_clim.nc')
fidclim = nc.Dataset('CCMP_1990-2010_clim.nc')

# get clim values
climtime  = fidclim.variables['time']
climlat   = fidclim.variables['latitude']
climlon   = fidclim.variables['longitude']
uclim = fidclim.variables['uwnd']
vclim = fidclim.variables['vwnd']

spval = vclim._Fillvalue

# change lon from 0 to 360 -> 0 to -180
lon_shift = climlon[:].copy()
lon_shift[climlon[:]>180] = climlon[climlon[:]>180]-2.*180

dir = '/glade/work/edrenkar/external_data/CCMP/'
ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
# LOOP OVER ALL YEARS
for year in range(2001,2002+1):
    print year
    # LEAP YEAR EVALUATION
    ineq = year%4
    if ineq == 0:
       ndays[1] = 29
    else:
       ndays[1] = 28

    # LOOP OVER ALL MONTHS
    for mon in range(12):

        # LOOP OVER ALL DAYS
        for day in range(ndays[mon]):
            ncfil= dir + str(year)+ '/CCMP_Wind_Analysis_'+ str(year) + str(mon+1).zfill(2)+str(day+1).zfill(2)+'_V02.0_L3.0_RSS.nc.sub'
            #print ncfil
            fid = nc.Dataset(ncfil)

            # DAILY
            # get 6-hrly time values; create daily means
            #tavg = np.mean(fid.variables['time'][:])
            #print tavg
            # get 6-hrly wind values; create daily means
            #uavg = np.mean(fid.variables['uwnd'][:],axis=0)
            #vavg = np.mean(fid.variables['vwnd'][:],axis=0)
            # generate wind field anomalies
            #vanom = vavg - vclim[mon,:,:].squeeze() 
            #uanom = uavg - uclim[mon,:,:].squeeze()

            # 6 HOURLY
            tim_hr = fid.variables['time'][:]
            uwind = fid.variables['uwnd'][:]
            vwind = fid.variables['vwnd'][:] 
            # generate wind field anomalies
            vanom = vwind - np.tile(vclim[mon,:,:].squeeze(),(4,1,1)) 
            uanom = uwind - np.tile(uclim[mon,:,:].squeeze(),(4,1,1))

            # create new anomaly file
            #outfil = dir + str(year)+ '/CCMP_Wind_Analysis_'+ str(year) + str(mon+1).zfill(2)+str(day+1).zfill(2) + '_daily_anom.nc'
            outfil = dir + str(year)+ '/CCMP_Wind_Analysis_'+ str(year) + str(mon+1).zfill(2)+str(day+1).zfill(2) + '_6hrly_anom.nc'
            newnc = nc.Dataset(outfil, 'w', format='NETCDF3_CLASSIC')

            # add global attributes
            for atname in fid.ncattrs():
                setattr(newnc,atname,getattr(fid,atname))

            fid.close()

            # add dimensions
            newnc.createDimension('time', None)
            newnc.createDimension('lat', len(climlat[:]))
            newnc.createDimension('lon', len(climlon[:]))

            time = newnc.createVariable('time', 'f8', ('time',))
            time.setncattr('units','Days')
            time.setncattr('calendar','NOLEAP')
            time[:4] = tim_hr/24.- (dt.datetime(year,1,1)-dt.datetime(1987,1,1)).days      

            lat = newnc.createVariable('lat', 'f8', ('lat',)) 
            lat.setncattr('units','degrees_north')
            lat.setncattr('long_name','latitude')
            lat[:] = climlat[:]
 
            lon = newnc.createVariable('lon', 'f8', ('lon',))
            lon.setncattr('units','degrees_east')
            lon.setncattr('long_name','longitude')
            lon[:] = lon_shift

            u_anom = newnc.createVariable( 'u_anom', 'f', ('time', 'lat', 'lon',))
            for atname in uclim.ncattrs():
                setattr(u_anom,atname,getattr(uclim,atname))
            # update certain attributes
            u_anom.standard_name = "eastward_wind_anomaly"
            u_anom.long_name = "u-wind vector component anomaly relative to 1990-2010 climatology at 10 meters"
            u_anom.valid_min = np.min(uanom)
            u_anom.valid_max = np.max(uanom)
            u_anom[:4] = uanom
            #u_anom[:1] = uanom

            v_anom = newnc.createVariable( 'v_anom', 'f', ('time', 'lat', 'lon',))
            for atname in vclim.ncattrs():
                setattr(v_anom,atname,getattr(vclim,atname)) 
            # update certain attributes
            v_anom.standard_name = "northward_wind_anomaly"
            v_anom.long_name = "v-wind vector component anomaly relative to 1990-2010 climatology at 10 meters"
            v_anom.valid_min = np.min(vanom)
            v_anom.valid_max = np.max(vanom)
            v_anom[:4] = vanom
            #v_anom[:1] = vanom 
  
            newnc.close()

fidclim.close()
