import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import ESMF

# MacBook
fid_ERAi = nc.Dataset('/Users/elizabethdrenkard/TOOLS/CCS_scripts/CCMP/ERAinterim_.125_lsm.nc')

elat = fid_ERAi.variables['lat'][:]
elon = fid_ERAi.variables['lon'][:]

elon = fid_ERAi.variables['lon'][:]
elon[elon>180]=elon[elon>180]-360
elat = fid_ERAi.variables['lat'][:]

eLSM = fid_ERAi.variables['LSM'][:]

Xi,Yi = np.meshgrid(elon,elat)
sourcegrid = ESMF.Grid(np.array(Yi.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

# CCMP file
yr = 'MAY01-APR02'
# MACBOOK
CCMP_fil = '/Users/elizabethdrenkard/external_data/CCMP/CCMP_'+yr+'_daily_anom.nc'  

fidCCMP = nc.Dataset(CCMP_fil)

clon = fidCCMP.variables['longitude'][:]
clon[clon>180]=clon[clon>180]-360
clat = fidCCMP.variables['latitude'][:]

nx = len(clon)
ny = len(clat)

Xn, Yn = np.meshgrid(clon,clat)
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

source_LSM = ESMF.Field(sourcegrid, name = 'ERAi_LSM')
dest_LSM = ESMF.Field(destgrid, name = 'CCMP_LSM')

source_LSM.data[...] = eLSM.squeeze()

regrid_LSM = ESMF.Regrid(source_LSM, dest_LSM, regrid_method = ESMF.RegridMethod.BILINEAR,  
                     unmapped_action = ESMF.UnmappedAction.IGNORE)

dest_LSM = regrid_LSM(source_LSM, dest_LSM) 

tmp_LSM = dest_LSM.data

# convert near-zero values to zero
tmp_LSM[tmp_LSM<.5] = 0

# reverse land/ocean mask values
new_LSM = np.zeros(tmp_LSM.shape)
new_LSM[tmp_LSM==0]=1

# Save new LSM file
ncfil_new = 'CCMP_LSM.nc'
fid2 = nc.Dataset(ncfil_new,'w')

fid2.createDimension('latitude', ny)
fid2.createDimension('longitude', nx)

fid2.createVariable('latitude','f8',('latitude'))
fid2.variables['latitude'].long_name = fid_ERAi.variables['lat'].long_name
fid2.variables['latitude'].units = fid_ERAi.variables['lat'].units
fid2.variables['latitude'][:]=clat

fid2.createVariable('longitude','f8',('longitude'))
fid2.variables['longitude'].long_name = fid_ERAi.variables['lon'].long_name
fid2.variables['longitude'].units = fid_ERAi.variables['lon'].units
fid2.variables['longitude'][:]=clon
    
fid2.createVariable('LSM','f8',('latitude','longitude'),fill_value = np.float(1.0e15))
fid2.variables['LSM'].long_name = fid_ERAi.variables['LSM'].long_name
fid2.variables['LSM'].code = fid_ERAi.variables['LSM'].code
fid2.variables['LSM'].table = fid_ERAi.variables['LSM'].table
u_txt = "Land-Sea Mask for CCMP derived from ERAinterim .125 degrees LSM"
fid2.variables['LSM'].details = u_txt
fid2.variables['LSM'][:]= new_LSM

fid2.close()

