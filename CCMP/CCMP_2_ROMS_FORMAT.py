import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# CCMP files
yr = 'MAY01-APR02'
# cheyenne
fid_ccmp01 = nc.Dataset('/glade/work/edrenkar/external_data/CCMP/CCMP_MAY-DEC01.nc')
fid_ccmp02 = nc.Dataset('/glade/work/edrenkar/external_data/CCMP/CCMP_JAN-APR02.nc')

U_01 = fid_ccmp01.variables['uwnd'][:]
U_02 = fid_ccmp02.variables['uwnd'][:]

V_01 = fid_ccmp01.variables['vwnd'][:]
V_02 = fid_ccmp02.variables['vwnd'][:]

# OLD ERAi FILES for metadata
fidU_OLD = nc.Dataset('/glade/work/edrenkar/external_data/CCMP/ERAi_CCMPanom_MAY01-APR02_Uwind.nc')
fidV_OLD = nc.Dataset('/glade/work/edrenkar/external_data/CCMP/ERAi_CCMPanom_MAY01-APR02_Vwind.nc')

# GEOGRAPHY
ccmp_lon = fid_ccmp01.variables['longitude'][:]
ccmp_lon[ccmp_lon>180]=ccmp_lon[ccmp_lon>180]-360
ccmp_lat = fid_ccmp01.variables['latitude'][:]

nx = len(ccmp_lon)
ny = len(ccmp_lat)

# Allocate output variables
Uout = np.append(U_02,U_01,axis=0)
Vout = np.append(V_02,V_01,axis=0) 

# Save new wind files
ncfilU = 'CCMP_' + yr + '_Uwind.nc'
fidU = nc.Dataset(ncfilU,'w')

fidU.createDimension('time', None)
fidU.createDimension('lat', ny)
fidU.createDimension('lon', nx)

fidU.createVariable('time', 'f8', ('time'))
fidU.variables['time'].units = fidU_OLD.variables['time'].units
fidU.variables['time'].cycle_length = fidU_OLD.variables['time'].cycle_length
fidU.variables['time'][:] = np.arange(0,365,.25)

fidU.createVariable('lat','f8',('lat'))
fidU.variables['lat'].long_name = fidU_OLD.variables['lat'].long_name
fidU.variables['lat'].units = fidU_OLD.variables['lat'].units
fidU.variables['lat'][:] = ccmp_lat

fidU.createVariable('lon','f8',('lon'))
fidU.variables['lon'].long_name = fidU_OLD.variables['lon'].long_name
fidU.variables['lon'].units = fidU_OLD.variables['lon'].units
fidU.variables['lon'][:] = ccmp_lon
    
fidU.createVariable('Uwind','f8',('time','lat','lon'),fill_value = np.float(1.0e15))
fidU.variables['Uwind'].long_name = fidU_OLD.variables['Uwind'].long_name
fidU.variables['Uwind'].units = fidU_OLD.variables['Uwind'].units
fidU.variables['Uwind'].coordinates = fidU_OLD.variables['Uwind'].coordinates
fidU.variables['Uwind'].time = fidU_OLD.variables['Uwind'].time
u_txt = "CCMP 6-hrly zonal wind component, " + yr
fidU.variables['Uwind'].details = u_txt
fidU.variables['Uwind'][:] = Uout

fidU.close()

ncfilV = 'CCMP_' + yr + '_Vwind.nc'
fidV = nc.Dataset(ncfilV,'w')
    
fidV.createDimension('time', None)
fidV.createDimension('lat',ny) 
fidV.createDimension('lon',nx) 
    
fidV.createVariable('time', 'f8', ('time'))
fidV.variables['time'].units = fidV_OLD.variables['time'].units
fidV.variables['time'].cycle_length = fidV_OLD.variables['time'].cycle_length
fidV.variables['time'][:] = np.arange(0,365,.25)

fidV.createVariable('lat','f8',('lat'))
fidV.variables['lat'].long_name = fidV_OLD.variables['lat'].long_name
fidV.variables['lat'].units = fidV_OLD.variables['lat'].units
fidV.variables['lat'][:]= ccmp_lat
    
fidV.createVariable('lon','f8',('lon'))
fidV.variables['lon'].long_name = fidV_OLD.variables['lon'].long_name
fidV.variables['lon'].units = fidV_OLD.variables['lon'].units
fidV.variables['lon'][:]= ccmp_lon
    
fidV.createVariable('Vwind','f8',('time','lat','lon'),fill_value = np.float(1.0e15))
fidV.variables['Vwind'].long_name = fidV_OLD.variables['Vwind'].long_name
fidV.variables['Vwind'].units = fidV_OLD.variables['Vwind'].units
fidV.variables['Vwind'].coordinates = fidV_OLD.variables['Vwind'].coordinates
fidV.variables['Vwind'].time = fidV_OLD.variables['Vwind'].time
v_txt = "CCMP 6-hrly meridional wind component, " + yr

fidV.variables['Vwind'].details = v_txt 
fidV.variables['Vwind'][:]= Vout

fidV.close()

fidV_OLD.close()
fidU_OLD.close()
