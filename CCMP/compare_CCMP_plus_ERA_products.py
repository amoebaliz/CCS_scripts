import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import ESMF

# Checking how drowned CCMP + ERAi compares with original, 
# non-drowned CCMP + ERAi

# NOTE: new version interpolates to the ERAi grid. 
#       old version compares to the CCMP grid... 

# METHOD: Interpolate NEW composite file  back to CCMP domain 
#         and subtract/visualize fields

# 0 for uwind   
# 1 for vwind
 
nvar = 1 

old_fils = ['ERAi_CCMPanom_MAY01-APR02_Uwind.nc.old','ERAi_CCMPanom_MAY01-APR02_Vwind.nc.old']
new_fils = ['ERAi_CCMPanom_MAY01-APR02_Uwind.nc','ERAi_CCMPanom_MAY01-APR02_Vwind.nc']

fvars = ['Uwind','Vwind'] 

fid_0 = nc.Dataset(old_fils[nvar])
fid_1 = nc.Dataset(new_fils[nvar])

# ORIGINAL FILES: on the sub CCMP domain
var_0 = fid_0.variables[fvars[nvar]][:]
lat_0 = fid_0.variables['lat'][:]
lon_0 = fid_0.variables['lon'][:]

Xn,Yn = np.meshgrid(lon_0,lat_0)
destgrid = ESMF.Grid(np.array(Yn.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

# NEW FILES: on the ERAi global domain
var_1 = fid_1.variables[fvars[nvar]][:]
lat_1 = fid_1.variables['lat'][:]
lon_1 = fid_1.variables['lon'][:]

Xi,Yi = np.meshgrid(lon_1,lat_1)
sourcegrid = ESMF.Grid(np.array(Yi.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

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

sourcefield = ESMF.Field(sourcegrid, name = 'ERAi+CCMPanom_ERAgrid')
destfield = ESMF.Field(destgrid, name = 'ERAi+CCMPanom_CCMPgrid')

for nt in range(lon_0.shape[0]):

    sourcefield.data[...] = var_1[nt,:].squeeze() 
    # REGULAR GRID BILINIEAR INTERPOLATION
    regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.BILINEAR,
                   unmapped_action = ESMF.UnmappedAction.IGNORE)

    destfield = regrid(sourcefield, destfield)

    out = var_0[nt+245,:].squeeze() - destfield.data
    #print np.sum(out)
    plt.pcolor(out) 
    plt.colorbar()
    plt.show()

  




