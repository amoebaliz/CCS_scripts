import numpy as np
import netCDF4 as nc
import ESMF


# ROMS Grid Information
grdfile = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid_roms = nc.Dataset(grdfile)
roms_lat = fid_roms.variables['lat_rho'][:]
roms_lon = fid_roms.variables['lon_rho'][:]
roms_mask = fid_roms.variables['mask_rho'][:]
roms_mask = -1*roms_mask+1

ny = roms_lat.shape[0]
nx = roms_lon.shape[1]

# FUT FILES: Reshape to CCS grid
fid1 = nc.Dataset('diatChl_may_delta.nc')
fid2 = nc.Dataset('diazChl_may_delta.nc')
fid3 = nc.Dataset('spChl_may_delta.nc')

# get lat lon
lon = fid1.variables['TLONG'][:]
lat = fid1.variables['TLAT'][:]
chl1 = fid1.variables['diatChl'][:]
chl2 = fid2.variables['diazChl'][:]
chl3 = fid3.variables['spChl'][:]

fut_chl_MAM = (chl1+chl2+chl3).squeeze()

print  fut_chl_MAM.shape


# FOR 1D lat/lon grids
sourcegrid = ESMF.Grid(np.array(lat.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

destgrid = ESMF.Grid(np.array(roms_lon.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

## POINTERS
source_lon = sourcegrid.get_coords(0)
source_lat = sourcegrid.get_coords(1)
dest_lon = destgrid.get_coords(0)
dest_lat = destgrid.get_coords(1)

## FILLS
source_lon[...] = lon 
source_lat[...] = lat
dest_lon[...] = roms_lon
dest_lat[...] = roms_lat

sourcefield = ESMF.Field(sourcegrid, name = 'CESM_chl')
destfield = ESMF.Field(destgrid, name = 'CCS_chl')

# REGULAR GRID BILINIEAR INTERPOLATION
sourcefield.data[...] = fut_chl_MAM.squeeze()
regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.BILINEAR,unmapped_action = ESMF.UnmappedAction.IGNORE)
destfield = regrid(sourcefield, destfield)

out_var = destfield.data
out_var[roms_mask.astype(np.bool)]=np.ma.masked
out_var[out_var>100]=np.ma.masked


# FOR METADATA ONLY
wifs_nc = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/SeaWIFS/SeaWiFS_CCS_his_Clim_MAM.nc'
fid4 = nc.Dataset(wifs_nc)


# Save new files
ncfile = 'CESMfut017_CCS_Clim_may.nc.full_rho'
fid2 = nc.Dataset(ncfile,'w')

fid2.createDimension('lat', ny)
fid2.createDimension('lon', nx)

fid2.createVariable('lat','f8',('lat','lon'))
fid2.variables['lat'].long_name = 'latitude'
fid2.variables['lat'][:]=roms_lat

fid2.createVariable('lon','f8',('lat','lon'))
fid2.variables['lon'].long_name = 'longitude'
fid2.variables['lon'][:]=roms_lon

fid2.createVariable('Chla','f8',('lat','lon'),fill_value = np.float(1.0e15))
fid2.variables['Chla'].long_name = fid4.variables['Chla'].long_name
fid2.variables['Chla'].units = fid4.variables['Chla'].units
u_txt = "mass_chlorophyll_concentration_in_sea_water"
fid2.variables['Chla'].standard_name = u_txt
fid2.variables['Chla'][:]= out_var 

fid2.close()

