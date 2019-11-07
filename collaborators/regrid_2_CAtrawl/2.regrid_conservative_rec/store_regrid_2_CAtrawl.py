import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import ESMF
import shapefile as shp 
import matplotlib.pyplot as plt

def get_var(i):
    if roms:
       # CONSTRUCT FILE NAME STRINGS 
       nc1 = his_dir + 'CCS-LD.HCo02Y_clim_' + str(i+1).zfill(2) + '.nc'
       nc2 = fut_dir + 'CCS-LD.FCo017_clim_' + str(i+1).zfill(2) + '.nc' 

       # IF MONTHLY FILES, NO NEED TO SPECIFY TIME INDEX
       fid1 = nc.Dataset(nc1)
       var1 = fid1.variables[ocean_var][:,dep_lev,:].squeeze()
       fid2 = nc.Dataset(nc2)
       var2 = fid2.variables[ocean_var][:,dep_lev,:].squeeze()
       var_dif = (var2-var1)[1:-1,1:-1] 

    else:
       fid = nc.Dataset(ncfile)
       var_dif = fid.variables[ocean_var][i,dep_lev,1:-1,1:-1].squeeze() 

    return var_dif

def interp_to_caltrawl(ns,destfield):
    print 'CELL # ', ns
    clim_vals = []
    for mon in range(nmon):
 
        var = get_var(mon)
        sourcefield.data[...] = var.squeeze()
        # INTERPOLATE TO SINGLE BOX 
        destfield = regrid(sourcefield, destfield)

        # EXTRACT CENTER VALUE
        c_val = destfield.data[...][0][0] 

        # IF NO OCEAN CELL OVERLAP, SET TO NAN
        if np.sum(srcfracfield.data) == 0:
           c_val=np.nan

        #print ns, mon_list[mon], c_val
        clim_vals.append(c_val)

    # Write climatology to shapefile  
    w.records.append(sf.records()[ns][:1] + clim_vals)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CALIFORNIA TRAWL INFORMATION
shp_fil = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/collaborators/Farrah_shp/caltrawl_GCS'
sf = shp.Reader(shp_fil)
nshapes = len(sf.shapes())
roms=0
# NEW SHAPEFILE
#w = shp.Writer(sf.shapeType)

# WRITE FIELDS TO NEW SHAPEFILE
mon_list =['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
nmon = 12
var_str = '_dSST' 
w.fields = [sf.fields[n] for n in range(2)]
for nm in range(nmon):
    w.fields.append([mon_list[nm]+var_str] + ['N', 19, 11])

# COPY SHAPEFILE SHAPES
w._shapes.extend(sf.shapes())

# ROMS GRID INFO FOR INTERPOLATION
grdfile = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid = nc.Dataset(grdfile)
if roms:
   mask_rho = fid.variables['mask_rho'][1:-1,1:-1]
   plat = fid.variables['lat_psi'][:]
   plon = fid.variables['lon_psi'][:]

else:
   # CESM GRID INFO FOR INTERPOLATION
   grd_fil = '/Users/elizabethdrenkard/Documents/Conferences/2018/ECCWO/ECCWO_FILES/LENS_grid.nc'
   fid = nc.Dataset(grd_fil)
   mask_rho = fid.variables['lsm'][:][1:-1,1:-1]
   print mask_rho.shape, mask_rho
   lats = fid.variables['lat'][:]
   print lats.shape
   plat = (lats[:-1,:-1] + lats[1:,:-1] + lats [1:,1:] + lats[:-1,1:])/4. 
   lons = fid.variables['lon'][:]
   lons[lons>180]=lons[lons>180]-360
   plon = (lons[:-1,:-1] + lons[1:,:-1] + lons [1:,1:] + lons[:-1,1:])/4.

# PREPARE SOURCE GRID OBJECT
sourcegrid = ESMF.Grid(np.array(mask_rho.shape), staggerloc = ESMF.StaggerLoc.CORNER, \
                       coord_sys = ESMF.CoordSys.SPH_DEG)

source_lon = sourcegrid.get_coords(0, staggerloc=ESMF.StaggerLoc.CORNER)
source_lat = sourcegrid.get_coords(1, staggerloc=ESMF.StaggerLoc.CORNER)

source_lon[...] = plon
source_lat[...] = plat

# MASKING LAND POINTS
sourcegrid.add_item(ESMF.GridItem.MASK,[ESMF.StaggerLoc.CENTER])
grid_mask = sourcegrid.get_item(ESMF.GridItem.MASK)
grid_mask[...] = mask_rho.astype(np.int32)

sourcefield = ESMF.Field(sourcegrid, name = 'ROMS_Delta')
srcfracfield = ESMF.Field(sourcegrid, 'srcfracfield')

# PREPARE DESTINATION GRID OBJECT
destgrid = ESMF.Grid(np.array((1,1)), staggerloc = ESMF.StaggerLoc.CORNER, \
                     coord_sys = ESMF.CoordSys.SPH_DEG)
dest_lon = destgrid.get_coords(0,staggerloc=ESMF.StaggerLoc.CORNER)
dest_lat = destgrid.get_coords(1,staggerloc=ESMF.StaggerLoc.CORNER)
if roms:
   # DEFINE ROMS OUTPUT INFORMATION
   model = 'ROMS'
   his_dir='/Volumes/Abalone/CCS/his2/climatologies/'
   fut_dir='/Volumes/Abalone/CCS/fut_017/climatologies/'
   ocean_var='temp'
   dep_lev=0 # Surface=49, Bottom=0
else:
   # DEFINE CESM OUTPUT INFORMATION
   model = 'CESM'
   ncfile = '/Volumes/Abalone/CCS/CESM/climatologies/017_TEMP_clim_delta.nc'
   ocean_var='TEMP'
   dep_lev=0 # Surface=0, Bottom=?

# LOOP OVER EACH SHAPE IN SHAPEFILE
w.records = []
for ns in range(nshapes):
    s = sf.shape(ns)
    bbox = [coord for coord in s.bbox]

    # EXTRACT CORNER COORDINATES
    shp_lons = np.array((bbox[0], bbox[2]))
    shp_lats = np.array((bbox[3], bbox[1]))

    Xn, Yn = np.meshgrid(shp_lons, shp_lats)

    # ASSIGN TO DESTINATION GRID OBJECT
    dest_lon[...] = Xn
    dest_lat[...] = Yn
    # ASSIGN DESTINATION OBJECT
    destfield = ESMF.Field(destgrid, name = 'CAtrawl_delta')

    # DEFINE INTERPOLATION FUNCTION
    regrid = ESMF.Regrid(sourcefield, destfield, 
                         regrid_method = ESMF.RegridMethod.CONSERVE,
                         src_mask_values=np.array([0], dtype=np.int32),
                         src_frac_field=srcfracfield,
                         norm_type=ESMF.NormType.FRACAREA, 
                         unmapped_action = ESMF.UnmappedAction.IGNORE)

    interp_to_caltrawl(ns,destfield)

# SAVE TO SHAPE FILE
new_shp_file = 'caltrawl_GCS_Climatological_' + model + var_str
w.save(new_shp_file)
