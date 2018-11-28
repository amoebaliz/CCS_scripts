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

# HISTORICAL FILES: Reshape to CCS grid
fid1 = nc.Dataset('SeaWiFS_Chla_03.nc')
fid2 = nc.Dataset('SeaWiFS_Chla_04.nc')
fid3 = nc.Dataset('SeaWiFS_Chla_05.nc')

# get lat lon
lon = fid1.variables['lon'][:]
lat = fid1.variables['lat'][:]
chl1 = fid1.variables['chlor_a'][:]
chl2 = fid2.variables['chlor_a'][:]
chl3 = fid3.variables['chlor_a'][:]


his_chl_MAM = (chl1+chl2+chl3)/3.

print  his_chl_MAM.shape


# FOR 1D lat/lon grids
Xi,Yi = np.meshgrid(lon,lat)
sourcegrid = ESMF.Grid(np.array(Xi.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

destgrid = ESMF.Grid(np.array(roms_lon.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

## POINTERS
source_lon = sourcegrid.get_coords(0)
source_lat = sourcegrid.get_coords(1)
dest_lon = destgrid.get_coords(0)
dest_lat = destgrid.get_coords(1)

## FILLS
source_lon[...] = Xi 
source_lat[...] = Yi
dest_lon[...] = roms_lon
dest_lat[...] = roms_lat

sourcefield = ESMF.Field(sourcegrid, name = 'seawifs_chl')
destfield = ESMF.Field(destgrid, name = 'SODA_chl')

# REGULAR GRID BILINIEAR INTERPOLATION
sourcefield.data[...] = his_chl_MAM.squeeze()
regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.BILINEAR,unmapped_action = ESMF.UnmappedAction.IGNORE)
destfield = regrid(sourcefield, destfield)

out_var = destfield.data
out_var[roms_mask.astype(np.bool)]=np.ma.masked
out_var[out_var<0]=0

# Save new files
ncfile = 'SeaWiFS_CCS_his_Clim_MAM.nc'
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
fid2.variables['Chla'].long_name = fid1.variables['chlor_a'].long_name
fid2.variables['Chla'].units = fid1.variables['chlor_a'].units
u_txt = "mass_concentration_chlorophyll_concentration_in_sea_water"
fid2.variables['Chla'].standard_name = u_txt
fid2.variables['Chla'][:]= out_var

fid2.close()

