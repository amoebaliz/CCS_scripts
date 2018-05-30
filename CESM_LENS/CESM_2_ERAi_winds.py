import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import ESMF




# ERAint + CCMP files
# MacBook file location
fidU_EC = nc.Dataset('ERAi_CCMPanom_MAY01-APR02_Uwind.nc')
fidV_EC = nc.Dataset('ERAi_CCMPanom_MAY01-APR02_Vwind.nc')

U_EC = fidU_EC.variables['Uwind'][:]
V_EC = fidV_EC.variables['Vwind'][:]

dlon = fidU_EC.variables['lon'][:]
dlat = fidU_EC.variables['lat'][:]

ny = len(dlat)
nx = len(dlon)

Xn,Yn = np.meshgrid(dlon,dlat)
destgrid = ESMF.Grid(np.array(Xn.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

yr = 'MAY01-APR02'

# CESM LENS DELTA
fidU_delta = nc.Dataset('/Users/elizabethdrenkard/external_data/CESM/drowned/drowned_LENS_017_u10_delta.nc')
fidV_delta = nc.Dataset('/Users/elizabethdrenkard/external_data/CESM/drowned/drowned_LENS_017_v10_delta.nc')

Udelta = fidU_delta.variables['Uwind'][:]
Vdelta = fidV_delta.variables['Vwind'][:]

lon = fidU_delta.variables['lon'][:]
lon[lon>180]=lon[lon>180]-360
lat = fidU_delta.variables['lat'][:]

Xi, Yi = np.meshgrid(lon,lat)
sourcegrid = ESMF.Grid(np.array(Xi.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

## POINTERS
source_lon = sourcegrid.get_coords(0)
source_lat = sourcegrid.get_coords(1)
dest_lon = destgrid.get_coords(0)
dest_lat = destgrid.get_coords(1)
## FILLS
source_lon[...] = Xi
source_lat[...] = Yi
dest_lon[...] = Xn
dest_lat[...] = Yn

sourcefieldU = ESMF.Field(sourcegrid, name = 'DESM_Uwind')
sourcefieldV = ESMF.Field(sourcegrid, name = 'CESM_DELTA_Vwind')

destfieldU = ESMF.Field(destgrid, name = 'ERA_CCMP_Anom')
destfieldV = ESMF.Field(destgrid, name = 'ERA_CCMP_Anom')
# Allocate output variables
Uout = np.zeros(U_EC.shape)
Vout = np.zeros(V_EC.shape)

## ITERATE OVER ALL DAYS TO ADD MONTHLY CLIM   
ndays = np.array((31,28,31,30,31,30,31,31,30,31,30,31))
n = 0

for nmon in range(12):
    # REGULAR GRID BILINIEAR INTERPOLATION
    sourcefieldU.data[...] = Udelta[nmon,:].squeeze()
    sourcefieldV.data[...] = Vdelta[nmon,:].squeeze()
    regridU = ESMF.Regrid(sourcefieldU, destfieldU, regrid_method = ESMF.RegridMethod.BILINEAR,  
                     unmapped_action = ESMF.UnmappedAction.IGNORE)
    regridV = ESMF.Regrid(sourcefieldV, destfieldV, regrid_method = ESMF.RegridMethod.BILINEAR,
                     unmapped_action = ESMF.UnmappedAction.IGNORE)
    destfieldU = regridU(sourcefieldU, destfieldU) 
    destfieldV = regridV(sourcefieldV, destfieldV)    
 
    for nt in range(ndays[nmon]):
        # daily index in CCMP file

        Uout[n,:] = U_EC[n,:].squeeze() + destfieldU.data
        Vout[n,:] = V_EC[n,:].squeeze() + destfieldV.data
        n+=1

# Save new wind files
ncfilU = 'ERAi_CCMPanom_CESM_017_delta_Uwind.nc'
fid2 = nc.Dataset(ncfilU,'w')

fid2.createDimension('time', None)
fid2.createDimension('lat', ny)
fid2.createDimension('lon', nx)

fid2.createVariable('time', 'f8', ('time'))
fid2.variables['time'].units = fidU_EC.variables['time'].units
fid2.variables['time'].cycle_length = fidU_EC.variables['time'].cycle_length
fid2.variables['time'][:] = np.arange(1,365+1)

fid2.createVariable('lat','f8',('lat'))
fid2.variables['lat'].long_name = fidU_EC.variables['lat'].long_name
fid2.variables['lat'].units = fidU_EC.variables['lat'].units
fid2.variables['lat'][:]=dlat

fid2.createVariable('lon','f8',('lon'))
fid2.variables['lon'].long_name = fidU_EC.variables['lon'].long_name
fid2.variables['lon'].units = fidU_EC.variables['lon'].units
fid2.variables['lon'][:]=dlon
    
fid2.createVariable('Uwind','f8',('time','lat','lon'),fill_value = np.float(1.0e15))
fid2.variables['Uwind'].long_name = fidU_EC.variables['Uwind'].long_name
fid2.variables['Uwind'].units = fidU_EC.variables['Uwind'].units
fid2.variables['Uwind'].coordinates = fidU_EC.variables['Uwind'].coordinates
fid2.variables['Uwind'].time = fidU_EC.variables['Uwind'].time
u_txt = "ERAinterim climatology (1981-2010) + CCMP " + yr + " Anomaly (relative to 1990-2010) + CESM LENS 017 DELTA"
fid2.variables['Uwind'].details = u_txt
fid2.variables['Uwind'][:]=Uout

fid2.close()

ncfilV = 'ERAi_CCMPanom_CESM_017_delta_Vwind.nc'
fid2 = nc.Dataset(ncfilV,'w')
    
fid2.createDimension('time', None)
fid2.createDimension('lat',ny) 
fid2.createDimension('lon',nx) 
    
fid2.createVariable('time', 'f8', ('time'))
fid2.variables['time'].units = fidV_EC.variables['time'].units
fid2.variables['time'].cycle_length = fidV_EC.variables['time'].cycle_length
fid2.variables['time'][:] = np.arange(1,365+1)    

fid2.createVariable('lat','f8',('lat'))
fid2.variables['lat'].long_name = fidV_EC.variables['lat'].long_name
fid2.variables['lat'].units = fidV_EC.variables['lat'].units
fid2.variables['lat'][:]=dlat
    
fid2.createVariable('lon','f8',('lon'))
fid2.variables['lon'].long_name = fidV_EC.variables['lon'].long_name
fid2.variables['lon'].units = fidV_EC.variables['lon'].units
fid2.variables['lon'][:]=dlon
    
fid2.createVariable('Vwind','f8',('time','lat','lon'),fill_value = np.float(1.0e15))
fid2.variables['Vwind'].long_name = fidV_EC.variables['Vwind'].long_name
fid2.variables['Vwind'].units = fidV_EC.variables['Vwind'].units
fid2.variables['Vwind'].coordinates = fidV_EC.variables['Vwind'].coordinates
fid2.variables['Vwind'].time = fidV_EC.variables['Vwind'].time
v_txt = "ERAinterim climatology (1981-2010) + CCMP " + yr + " Anomaly (relative to 1990-2010) + CESM LENS 017 DELTA"

fid2.variables['Vwind'].details = v_txt 
fid2.variables['Vwind'][:]=Vout

fid2.close()

