import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import ESMF

# ERAint files
# yellowstone
#fidUclim = nc.Dataset('/glade2/scratch2/edrenkar/CCS-inputs/drowned_ERAi_u10_1981-2010_monthly_clim.nc')
#fidVclim = nc.Dataset('/glade2/scratch2/edrenkar/CCS-inputs/drowned_ERAi_v10_1981-2010_monthly_clim.nc')

# MacBook

fidUdelta = nc.Dataset('/Users/elizabethdrenkard/Desktop/016/drowned_LENS_016_u10_delta.nc')
fidVdelta = nc.Dataset('/Users/elizabethdrenkard/Desktop/016/drowned_LENS_016_v10_delta.nc')

Udelta = fidUdelta.variables['Uwind'][:]
Vdelta = fidVdelta.variables['Vwind'][:]

dlon = fidUdelta.variables['lon'][:]
dlon[dlon>180]=dlon[dlon>180]-360
dlat = fidUdelta.variables['lat'][:]

Xi,Yi = np.meshgrid(dlon,dlat)
sourcegrid = ESMF.Grid(np.array(Xi.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

# CCMP file
yr = 'MAY01-APR02'
# MACBOOK
ERA_CCMP_ufil = '/Users/elizabethdrenkard/Desktop/016/ERAi_CCMPanom_MAY01-APR02_Uwind.nc'  
ERA_CCMP_vfil = '/Users/elizabethdrenkard/Desktop/016/ERAi_CCMPanom_MAY01-APR02_Vwind.nc' 
# YELLOWSTONE
#CCMP_fil = '/glade/p/work/edrenkar/external_data/CCMP/CCMP_'+yr+'_daily_anom.nc'   
fidCCMPU = nc.Dataset(ERA_CCMP_ufil)
fidCCMPV = nc.Dataset(ERA_CCMP_vfil)

FinU = fidCCMPU.variables['Uwind'][:]
FinV = fidCCMPV.variables['Vwind'][:]

lon = fidCCMPU.variables['lon'][:]
lon[lon>180]=lon[lon>180]-360
lat = fidCCMPU.variables['lat'][:]

nx = len(lon)
ny = len(lat)

Xn, Yn = np.meshgrid(lon,lat)
destgrid = ESMF.Grid(np.array(Xn.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)
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

sourcefieldU = ESMF.Field(sourcegrid, name = 'CCMP_DELTA_Uwind')
sourcefieldV = ESMF.Field(sourcegrid, name = 'CCMP_DELTA_Vwind')

destfieldU = ESMF.Field(destgrid, name = 'ERA_CCMP_Anom')
destfieldV = ESMF.Field(destgrid, name = 'ERA_CCMP_Anom')
# Allocate output variables
Uout = np.zeros(FinU.shape)
Vout = np.zeros(FinV.shape)

## ITERATE OVER ALL DAYS TO ADD MONTHLY CLIM   
ndays = np.array((31,28,31,30,31,30,31,31,30,31,30,31))
n = 0
for nmon in range(12):
    # REGULAR GRID BILINIEAR INTERPOLATION
    #regCU =  griddata(clon,clat,Uclim[nmon,:].squeeze(),lon,lat)
    #regCV =  griddata(clon,clat,Vclim[nmon,:].squeeze(),lon,lat)
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

        Uout[n,:] = FinU[n,:].squeeze() + destfieldU.data
        Vout[n,:] = FinV[n,:].squeeze() + destfieldV.data
        n+=1
# Save new wind files
ncfilU = 'ERAi_CCMPanom_CESM_016_delta_Uwind.nc'
fid2 = nc.Dataset(ncfilU,'w')

fid2.createDimension('time', None)
fid2.createDimension('lat', ny)
fid2.createDimension('lon', nx)

fid2.createVariable('time', 'f8', ('time'))
fid2.variables['time'].units = fidCCMPU.variables['time'].units
fid2.variables['time'].cycle_length = fidCCMPU.variables['time'].cycle_length
fid2.variables['time'][:] = np.arange(1,365+1)

fid2.createVariable('lat','f8',('lat'))
fid2.variables['lat'].long_name = fidCCMPU.variables['lat'].long_name
fid2.variables['lat'].units = fidCCMPU.variables['lat'].units
fid2.variables['lat'][:]=lat

fid2.createVariable('lon','f8',('lon'))
fid2.variables['lon'].long_name = fidCCMPU.variables['lon'].long_name
fid2.variables['lon'].units = fidCCMPU.variables['lon'].units
fid2.variables['lon'][:]=lon
    
fid2.createVariable('Uwind','f8',('time','lat','lon'),fill_value = np.float(1.0e15))
fid2.variables['Uwind'].long_name = fidCCMPU.variables['Uwind'].long_name
fid2.variables['Uwind'].units = fidCCMPU.variables['Uwind'].units
fid2.variables['Uwind'].coordinates = fidCCMPU.variables['Uwind'].coordinates
fid2.variables['Uwind'].time = fidCCMPU.variables['Uwind'].time
u_txt = "ERAinterim climatology (1981-2010) + CCMP " + yr + " Anomaly (relative to 1990-2010) + CESM LENS 016 DELTA"
fid2.variables['Uwind'].details = u_txt
fid2.variables['Uwind'][:]=Uout

fid2.close()

ncfilV = 'ERAi_CCMPanom_CESM_016_delta_Vwind.nc'
fid2 = nc.Dataset(ncfilV,'w')
    
fid2.createDimension('time', None)
fid2.createDimension('lat',ny) 
fid2.createDimension('lon',nx) 
    
fid2.createVariable('time', 'f8', ('time'))
fid2.variables['time'].units = fidCCMPV.variables['time'].units
fid2.variables['time'].cycle_length = fidCCMPV.variables['time'].cycle_length
fid2.variables['time'][:] = np.arange(1,365+1)    

fid2.createVariable('lat','f8',('lat'))
fid2.variables['lat'].long_name = fidCCMPV.variables['lat'].long_name
fid2.variables['lat'].units = fidCCMPV.variables['lat'].units
fid2.variables['lat'][:]=lat
    
fid2.createVariable('lon','f8',('lon'))
fid2.variables['lon'].long_name = fidCCMPV.variables['lon'].long_name
fid2.variables['lon'].units = fidCCMPV.variables['lon'].units
fid2.variables['lon'][:]=lon
    
fid2.createVariable('Vwind','f8',('time','lat','lon'),fill_value = np.float(1.0e15))
fid2.variables['Vwind'].long_name = fidCCMPV.variables['Vwind'].long_name
fid2.variables['Vwind'].units = fidCCMPV.variables['Vwind'].units
fid2.variables['Vwind'].coordinates = fidCCMPV.variables['Vwind'].coordinates
fid2.variables['Vwind'].time = fidCCMPV.variables['Vwind'].time
v_txt = "ERAinterim climatology (1981-2010) + CCMP " + yr + " Anomaly (relative to 1990-2010) + CESM LENS 016 DELTA"

fid2.variables['Vwind'].details = v_txt 
fid2.variables['Vwind'][:]=Vout

fid2.close()

fidUdelta.close()
fidVdelta.close()
