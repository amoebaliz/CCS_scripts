import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import shapefile as shp
import ESMF

# SOURCE GRID FILE 

# PREPARE SOURCE GRID OBJECT

# CALIFORNIA TRAWL INFORMATION
shp_fil = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/collaborators/Farrah_shp/caltrawl_GCS'
sf = shp.Reader(shp_fil)
nshapes = len(sf.shapes())

# PREPARE DESTINATION GRID OBJECT
destgrid = ESMF.Grid(np.array(3,3), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)
dest_lon = destgrid.get_coords(0)
dest_lat = destgrid.get_coords(1)

# LOOP OVER EACH SHAPE
for ns in range(nshapes):
    s = sf.shape(ns)
    bbox = [coord for coord in s.bbox]

    c_lon = np.mean((bbox[0], bbox[2]))
    c_lat = np.mean((bbox[3], bbox[1]))

    x = [2*bbox[0]-c_lon, c_lon, 2*bbox[2]-c_lon]
    y = [2*bbox[1]-c_lat, c_lat, 2*bbox[3]-c_lat]

    Xn, Yn = np.meshgrid(x, y)
    dest_lon[...] = Xn
    dest_lat[...] = Yn

    regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.BILINEAR,
             unmapped_action = ESMF.UnmappedAction.IGNORE)

    destfield = regrid(sourcefield, destfield)


