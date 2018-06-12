import pyroms
import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import matplotlib.pyplot as plt 
import matplotlib.dates as pltd

def get_salt():
    ncfile = '/Users/elizabethdrenkard/Desktop/HIS_5yr_SALT.nc'
    fid = nc.Dataset(ncfile)
    salt = fid.variables['salt'][:].squeeze()
    return salt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

grd = pyroms.grid.get_ROMS_grid('CCS')
#glat = grd.hgrid.lat_rho[1:-1,1:-1]
#glon = grd.hgrid.lon_rho[1:-1,1:-1]
#angs = grd.hgrid.angle_rho[1:-1,1:-1]
#x,y  = angs.shape 

istart = int(3)
iend   = int(302)
jstart = int(486)
jend   = int(326)

salt = get_salt()
#cs = get_vel()
#print cs.shape
#cs_transect, z, lon, lat = pyroms.tools.transect(cs, istart, iend, jstart, jend, grd,vert=False, Cpos='rho')

#x = np.tile(np.array(range(transect[:].shape[1])),(transect[:].shape[0],1))
#v = np.linspace(-1*col_val,col_val, 50, endpoint=True)

fig = plt.figure()
ax = fig.add_subplot(111,axisbg=[0.5,0.5,0.5])
v = np.linspace(-1*col_val,col_val, 50, endpoint=True)
cs = ax.contourf(x,z,transect,v,extend = 'both', cmap='bwr')
plt.colorbar(cs,ticks=[-.001, 0, 0.001])
plt.show()

