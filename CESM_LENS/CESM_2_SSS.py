import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import ESMF

# ERAint files
# yellowstone

# MacBook SSS file
sss_dir = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/SSS/'
sss_fil = sss_dir + 'sss_monthly_climatology_flooded.nc'
sss_fid = nc.Dataset(sss_fil)
sss_var = sss_fid.variables['SSS'][:]
slat = sss_fid.variables['lat'][:] 
slon = sss_fid.variables['lon'][:]
slon[slon>180]=slon[slon>180]-360

Xn, Yn = np.meshgrid(slon,slat)
destgrid = ESMF.Grid(np.array(Xn.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

# MacBook LENS salinity 
cmod = 16

lens_dir  = '/Users/elizabethdrenkard/Desktop/016/'
fil_delta = lens_dir + 'SSS_016_delta.nc' 
fid_delta = nc.Dataset(fil_delta) 
delta_var = fid_delta.variables['SALT'][:].squeeze()
dlon = fid_delta.variables['TLONG'][:]
dlon[dlon>180]=dlon[dlon>180]-360
dlat = sss_fid.variables['TLAT'][:]

#Xi,Yi = np.meshgrid(dlon,dlat)
#sourcegrid = ESMF.Grid(np.array(Xi.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)
sourcegrid = ESMF.Grid(np.array(dlon.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

# Date details
ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
time_vals = np.zeros(12)
dtot = 0 

for n in range(12):
    time_vals[n] = ndays[n]/2.0 + dtot
    dtot += ndays[n]

# Processing SSS LESNS variable

## POINTERS
source_lon = sourcegrid.get_coords(0)
source_lat = sourcegrid.get_coords(1)
dest_lon = destgrid.get_coords(0)
dest_lat = destgrid.get_coords(1)

## FILLS
source_lon[...] = dlon
source_lat[...] = dlat
dest_lon[...] = Xn
dest_lat[...] = Yn

sourcefield = ESMF.Field(sourcegrid, name = 'CESM_Delta')
destfield = ESMF.Field(destgrid, name = 'PHC3.0_Clim')

# Allocate output variables
var_out = np.zeros(sss_var.shape)

## ITERATE OVER ALL MONTHS 
for nmon in range(12):
    # REGULAR GRID BILINIEAR INTERPOLATION
    sourcefield.data[...] = lens_var[nmon,:].squeeze()
    regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.BILINEAR,  
             unmapped_action = ESMF.UnmappedAction.IGNORE)
    destfield = regrid(sourcefield, destfield) 
   
    var_out[nmon,:] = era_var[nmon,:].squeeze() + destfield.data
 
# Save new sss file
ncfile = 'SSS_CESM_' + str(cmod).zfill(3)+ '_delta.nc'
fid2 = nc.Dataset(ncfile,'w')

fid2.createDimension('sss_time', None)
fid2.createDimension('lat', ny)
fid2.createDimension('lon', nx)

fid2.createVariable('sss_time', 'f8', ('sss_time'))
fid2.variables['sss_time'].units = sss_fid.variables['sss_time'].units
fid2.variables['sss_time'].cycle_length = sss_fid.variables['sss_time'].cycle_length
fid2.variables['sss_time'][:] = time_vals

fid2.createVariable('lat','f8',('lat'))
fid2.variables['lat'].long_name = sss_fid.variables['lat'].long_name
fid2.variables['lat'].units = sss_fid.variables['lat'].units
fid2.variables['lat'][:]=lat

fid2.createVariable('lon','f8',('lon'))
fid2.variables['lon'].long_name = sss_fid.variables['lon'].long_name
fid2.variables['lon'].units = sss_fid.variables['lon'].units
fid2.variables['lon'][:]=lon
    
fid2.createVariable('SSS','f8',('sss_time','lat','lon'),fill_value = np.float(1.0e15))
fid2.variables['SSS'].long_name = sss_fid.variables[Var_nm[nv]].long_name
fid2.variables['SSS'].units = sss_fid.variables[Var_nm[nv]].units
fid2.variables['SSS'].coordinates = sss_fid.variables[Var_nm[nv]].coordinates
fid2.variables['SSS'].time = sss_fid.variables[Var_nm[nv]].time
u_txt = "PHC3.0 climatology + CESM Delta"
fid2.variables['SSS'].details = u_txt
fid2.variables['SSS'][:]=var_out

fid2.close()
sss_fid.close()
