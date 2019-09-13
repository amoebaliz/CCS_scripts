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
#fidUclim = nc.Dataset('/Users/elizabethdrenkard/external_data/ERAinterim/drowned/drowned_ERAi_u10_1981-2010_monthly_clim.nc')
#fidVclim = nc.Dataset('/Users/elizabethdrenkard/external_data/ERAinterim/drowned/drowned_ERAi_v10_1981-2010_monthly_clim.nc')

#Cheyenne
fidUclim = nc.Dataset('/glade/work/edrenkar/external_data/ERAinterim/drowned/6hry_or_daily/drowned_u10_ERAinterim_1981-2010_6hrly_clim.nc')
fidVclim = nc.Dataset('/glade/work/edrenkar/external_data/ERAinterim/drowned/6hry_or_daily/drowned_v10_ERAinterim_1981-2010_6hrly_clim.nc')

Uclim = fidUclim.variables['Uwind'][:]
Vclim = fidVclim.variables['Vwind'][:]

elon = fidUclim.variables['lon'][:]
elon[elon>180]=elon[elon>180]-360
elat = fidUclim.variables['lat'][:]

nx = len(elon)
ny = len(elat)

Xi,Yi = np.meshgrid(elon,elat)

# CCMP file
yr = 'MAY01-APR02'

# CHEYENNE
CCMP_u = '/glade/work/edrenkar/external_data/CCMP/drowned_CCMP_MAY01-APR02_Uwind.nc'
CCMP_v = '/glade/work/edrenkar/external_data/CCMP/drowned_CCMP_MAY01-APR02_Vwind.nc'

fidCCMPu = nc.Dataset(CCMP_u)
fidCCMPv = nc.Dataset(CCMP_v)

anomU = fidCCMPu.variables['Uwind'][:]
anomV = fidCCMPv.variables['Vwind'][:]

clon = fidCCMPu.variables['lon'][:]
clon[clon>180]=clon[clon>180]-360
clat = fidCCMPu.variables['lat'][:]

nx = len(clon)
ny = len(clat)

Xn, Yn = np.meshgrid(clon,clat)

# Source Grid = ERAinterim
# Destination Grid = CCMP
sourcegrid = ESMF.Grid(np.array(Xi.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)
destgrid = ESMF.Grid(np.array(Yn.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)
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

# CCMP to ERAi
sourcefieldU = ESMF.Field(sourcegrid, name = 'CCMP_Anom')
sourcefieldV = ESMF.Field(sourcegrid, name = 'CCMP_Anom')

destfieldU = ESMF.Field(destgrid, name = 'ERAi_CCMPanom_Uwind')
destfieldV = ESMF.Field(destgrid, name = 'ERAi_CCMPanom_Vwind')

# Allocate output variables
Uout = np.zeros((4*365,ny,nx))
Vout = np.zeros((4*365,ny,nx))

for nt in range(4*365):
    print nt
    sourcefieldU.data[...] = Uclim[nt,:].squeeze()
    sourcefieldV.data[...] = Vclim[nt,:].squeeze()

    # REGULAR GRID BILINIEAR INTERPOLATION
    regridU = ESMF.Regrid(sourcefieldU, destfieldU, regrid_method = ESMF.RegridMethod.BILINEAR,
              unmapped_action = ESMF.UnmappedAction.IGNORE)
    regridV = ESMF.Regrid(sourcefieldV, destfieldV, regrid_method = ESMF.RegridMethod.BILINEAR,
              unmapped_action = ESMF.UnmappedAction.IGNORE)

    destfieldU = regridU(sourcefieldU, destfieldU)
    destfieldV = regridV(sourcefieldV, destfieldV)

    Uout[nt,:] = anomU[nt,:].squeeze() + destfieldU.data
    Vout[nt,:] = anomV[nt,:].squeeze() + destfieldV.data

# Save new wind files
ncfilU = 'ERAi_CCMPanom_' + yr + '_Uwind.nc'
fid2 = nc.Dataset(ncfilU,'w')

fid2.createDimension('time', None)
fid2.createDimension('lat', ny)
fid2.createDimension('lon', nx)

time = fid2.createVariable('time', 'f8', ('time'))
fid2.variables['time'].units = fidUclim.variables['time'].units
fid2.variables['time'].cycle_length = fidUclim.variables['time'].cycle_length
fid2.variables['time'][:] = np.arange(.5,364.5+1)

fid2.createVariable('lat','f8',('lat'))
fid2.variables['lat'].long_name = fidUclim.variables['lat'].long_name
fid2.variables['lat'].units = fidUclim.variables['lat'].units
fid2.variables['lat'][:]=elat

fid2.createVariable('lon','f8',('lon'))
fid2.variables['lon'].long_name = fidUclim.variables['lon'].long_name
fid2.variables['lon'].units = fidUclim.variables['lon'].units
fid2.variables['lon'][:]=elon
    
fid2.createVariable('Uwind','f8',('time','lat','lon'),fill_value = np.float(1.0e15))
fid2.variables['Uwind'].long_name = fidUclim.variables['Uwind'].long_name
fid2.variables['Uwind'].units = fidUclim.variables['Uwind'].units
fid2.variables['Uwind'].coordinates = fidUclim.variables['Uwind'].coordinates
fid2.variables['Uwind'].time = fidUclim.variables['Uwind'].time
u_txt = "ERAinterim climatology (1981-2010) + CCMP " + yr + " Anomaly (relative to 1990-2010)"
fid2.variables['Uwind'].details = u_txt
fid2.variables['Uwind'][:]=Uout
print time
fid2.close()

ncfilV = 'ERAi_CCMPanom_' +yr+'_Vwind.nc'
fid2 = nc.Dataset(ncfilV,'w')
    
fid2.createDimension('time', None)
fid2.createDimension('lat',ny) 
fid2.createDimension('lon',nx) 
    
fid2.createVariable('time', 'f8', ('time'))
fid2.variables['time'].units = fidVclim.variables['time'].units
fid2.variables['time'].cycle_length = fidVclim.variables['time'].cycle_length
fid2.variables['time'][:] = np.arange(.5,364.5+1)    

fid2.createVariable('lat','f8',('lat'))
fid2.variables['lat'].long_name = fidVclim.variables['lat'].long_name
fid2.variables['lat'].units = fidVclim.variables['lat'].units
fid2.variables['lat'][:]=elat
    
fid2.createVariable('lon','f8',('lon'))
fid2.variables['lon'].long_name = fidVclim.variables['lon'].long_name
fid2.variables['lon'].units = fidVclim.variables['lon'].units
fid2.variables['lon'][:]=elon
    
fid2.createVariable('Vwind','f8',('time','lat','lon'),fill_value = np.float(1.0e15))
fid2.variables['Vwind'].long_name = fidVclim.variables['Vwind'].long_name
fid2.variables['Vwind'].units = fidVclim.variables['Vwind'].units
fid2.variables['Vwind'].coordinates = fidVclim.variables['Vwind'].coordinates
fid2.variables['Vwind'].time = fidVclim.variables['Vwind'].time
v_txt = "ERAinterim climatology (1981-2010) + CCMP " + yr + " Anomaly (relative to 1990-2010)"

fid2.variables['Vwind'].details = v_txt 
fid2.variables['Vwind'][:]=Vout

fid2.close()

fidVclim.close()
fidUclim.close()
