import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import ESMF
import shapefile as shp 

def get_var(i):

    # CONSTRUCT FILE NAME STRINGS 
    nc1 = his_dir + 'CCS-LD.HCo02Y_clim_' + str(i+1).zfill(2) + '.nc'
    nc2 = fut_dir + 'CCS-LD.FCo017_clim_' + str(i+1).zfill(2) + '.nc' 

    # IF MONTHLY FILES, NO NEED TO SPECIFY TIME INDEX
    fid1 = nc.Dataset(nc1)
    var1 = fid1.variables[ocean_var][:,dep_lev,:].squeeze()
    
    fid2 = nc.Dataset(nc2)
    var2 = fid2.variables[ocean_var][:,dep_lev,:].squeeze()
    var_dif = var2-var1 
    print np.min(var_dif),np.max(var_dif)
    return var_dif

def interp_to_caltrawl(var):
    im_all=[]
    sourcefield.data[...] = var.squeeze()
    # LOOP OVER EACH SHAPE IN SHAPEFILE
    for ns in range(nshapes):
    #for ns in range(20):
        s = sf.shape(ns)
        bbox = [coord for coord in s.bbox]
       
        # EXTRACT CORNER COORDINATES
        shp_lons = np.array((bbox[0], bbox[2]))
        shp_lats = np.array((bbox[3], bbox[1]))

        # DEFINE CENTRAL POINT
        c_lon = np.mean(shp_lons)
        c_lat = np.mean(shp_lats)
        # DEFINE SURROUNDING (FAKE) GRID POINTS
        x = [2*bbox[0]-c_lon, c_lon, 2*bbox[2]-c_lon]
        y = [2*bbox[1]-c_lat, c_lat, 2*bbox[3]-c_lat]

        Xn, Yn = np.meshgrid(x, y)
        dest_lon[...] = Xn
        dest_lat[...] = Yn
        
        # DEFINE INTERPOLATION FUNCTION
        destfield = ESMF.Field(destgrid, name = 'CAtrawl_delta')

        regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.BILINEAR,
                 unmapped_action = ESMF.UnmappedAction.IGNORE)
 
        # INTERPOLATE TO SINGLE BOX 
        destfield = regrid(sourcefield, destfield)

        # EXTRACT CENTER VALUE
        c_val = destfield.data[1,1] 
        #print ns, c_val
        if c_val>100:
           c_val=np.nan

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CALIFORNIA TRAWL INFORMATION
shp_fil = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/collaborators/Farrah_shp/caltrawl_GCS'
sf = shp.Reader(shp_fil)
nshapes = len(sf.shapes())

# ROMS Grid information
grdfile = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid = nc.Dataset(grdfile)
mask_rho = fid.variables['mask_rho'][:]
rlat = fid.variables['lat_rho'][:]
rlon = fid.variables['lon_rho'][:]
plat = fid.variables['lat_psi'][:]
plon = fid.variables['lon_psi'][:]
dep = fid.variables['h'][:]

# PREPARE SOURCE GRID OBJECT
sourcegrid = ESMF.Grid(np.array(mask_rho.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)
source_lon = sourcegrid.get_coords(0)
source_lat = sourcegrid.get_coords(1)
source_lon[...] = rlon
source_lat[...] = rlat
sourcefield = ESMF.Field(sourcegrid, name = 'ROMS_Delta')

# PREPARE DESTINATION GRID OBJECT
destgrid = ESMF.Grid(np.array((3,3)), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)
dest_lon = destgrid.get_coords(0)
dest_lat = destgrid.get_coords(1)
#destfield = ESMF.Field(destgrid, name = 'CAtrawl_delta')

# DEFINE ROMS OUTPUT INFORMATION
his_dir='/Volumes/Abalone/CCS/his2/climatologies/'
fut_dir='/Volumes/Abalone/CCS/fut_017/climatologies/'

ocean_var='temp'
dep_lev=49 # Surface=49, Bottom=0

for mon in range(12):
    var= get_var(i)
    interp_to_caltrawl(var)
    # SAVE TO SHAPE FILE
