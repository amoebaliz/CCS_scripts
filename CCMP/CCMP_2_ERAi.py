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
fidUclim = nc.Dataset('/Users/elizabethdrenkard/external_data/ERAinterim/drowned/drowned_ERAi_u10_1981-2010_monthly_clim.nc')
fidVclim = nc.Dataset('/Users/elizabethdrenkard/external_data/ERAinterim/drowned/drowned_ERAi_v10_1981-2010_monthly_clim.nc')

Uclim = fidUclim.variables['Uwind'][:]
Vclim = fidVclim.variables['Vwind'][:]

elon = fidUclim.variables['lon'][:]
elon[elon>180]=elon[elon>180]-360
elat = fidUclim.variables['lat'][:]

Xn,Yn = np.meshgrid(elon,elat)
destgrid = ESMF.Grid(np.array(Yn.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)


# CCMP file
mon0 = 5 # START MONTH (MAY) FOR CCMP
yr = 'MAY01-APR02'
# MACBOOK
CCMP_fil = '/Users/elizabethdrenkard/external_data/CCMP/CCMP_'+yr+'_daily_anom.nc'  
# YELLOWSTONE
#CCMP_fil = '/glade/p/work/edrenkar/external_data/CCMP/CCMP_'+yr+'_daily_anom.nc'   
fidCCMP = nc.Dataset(CCMP_fil)
anomU = fidCCMP.variables['u_anom'][:]
anomV = fidCCMP.variables['v_anom'][:]

clon = fidCCMP.variables['longitude'][:]
clon[clon>180]=clon[clon>180]-360
clat = fidCCMP.variables['latitude'][:]

nx = len(elon)
ny = len(elat)

Xi, Yi = np.meshgrid(clon,clat)
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

# CCMP to ERAi
sourcefieldU = ESMF.Field(sourcegrid, name = 'CCMP_Anom')
sourcefieldV = ESMF.Field(sourcegrid, name = 'CCMP_Anom')

destfieldU = ESMF.Field(destgrid, name = 'ERAi_CCMPanom_Uwind')
destfieldV = ESMF.Field(destgrid, name = 'ERAi_CCMPanom_Vwind')

# Allocate output variables
Uout = np.zeros((365,ny,nx))
Vout = np.zeros((365,ny,nx))

## ITERATE OVER ALL DAYS TO ADD MONTHLY CLIM   
## CCMP FILE (Starts with MAY, ends with APRIL)
ndays = np.array((31,30,31,31,30,31,30,31,31,28,31,30))
n2 = 0
for nmon in range(12):
 
    md = nmon-(mon0-1)
    nd = np.sum(ndays[:md])
    
    for nt in range(ndays[md]):        
        # daily index in CCMP file
        n = nd + nt
        print nmon, n, n2
        # add regridded clim to CCMP anomaly 
        sourcefieldU.data[...] = anomU[n,:].squeeze()
        sourcefieldV.data[...] = anomV[n,:].squeeze()

        # REGULAR GRID BILINIEAR INTERPOLATION
        regridU = ESMF.Regrid(sourcefieldU, destfieldU, regrid_method = ESMF.RegridMethod.BILINEAR,
                     unmapped_action = ESMF.UnmappedAction.IGNORE)
        regridV = ESMF.Regrid(sourcefieldV, destfieldV, regrid_method = ESMF.RegridMethod.BILINEAR,
                     unmapped_action = ESMF.UnmappedAction.IGNORE)

        destfieldU = regridU(sourcefieldU, destfieldU)
        destfieldV = regridV(sourcefieldV, destfieldV)

        Uout[n2,:] = Uclim[nmon,:].squeeze() + destfieldU.data
        Vout[n2,:] = Vclim[nmon,:].squeeze() + destfieldV.data

        n2+=1

# Save new wind files
ncfilU = 'ERAi_CCMPanom_' + yr + '_Uwind.nc'
fid2 = nc.Dataset(ncfilU,'w')

fid2.createDimension('time', None)
fid2.createDimension('lat', ny)
fid2.createDimension('lon', nx)

fid2.createVariable('time', 'f8', ('time'))
fid2.variables['time'].units = fidUclim.variables['time'].units
fid2.variables['time'].cycle_length = fidUclim.variables['time'].cycle_length
fid2.variables['time'][:] = np.arange(1,365+1)

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

fid2.close()

ncfilV = 'ERAi_CCMPanom_' +yr+'_Vwind.nc'
fid2 = nc.Dataset(ncfilV,'w')
    
fid2.createDimension('time', None)
fid2.createDimension('lat',ny) 
fid2.createDimension('lon',nx) 
    
fid2.createVariable('time', 'f8', ('time'))
fid2.variables['time'].units = fidVclim.variables['time'].units
fid2.variables['time'].cycle_length = fidVclim.variables['time'].cycle_length
fid2.variables['time'][:] = np.arange(1,365+1)    

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
