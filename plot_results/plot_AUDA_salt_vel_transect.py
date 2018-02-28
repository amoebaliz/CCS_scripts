import pyroms
import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import matplotlib.pyplot as plt 
import matplotlib.dates as pltd

def get_salt():
    ncfile = '/Users/elizabethdrenkard/Desktop/HIS_4yr_SALT.nc'
<<<<<<< HEAD
=======
    ncfile = '/glade/p/work/edrenkar/MODELS/CCS/RUNS/CCS-LD.HCo01Y/HIS_4yr_SALT.nc'
    ncfile = '/glade/p/work/edrenkar/MODELS/CCS/RUNS/4yr_clim_SST.nc' 
>>>>>>> 5a3e1fe4f1e70bd3874059569db971ddb881831c
    fid = nc.Dataset(ncfile)
    sst = fid.variables['salt'][:].squeeze()
    return sst[1:-1,1:-1]

def get_vel():
    ncfile = '/Users/elizabethdrenkard/Desktop/HIS_4yr_UV.nc'
    store = np.zeros((50,x+2,y+2))
    fid = nc.Dataset(ncfile)
    u_vel = np.squeeze(fid.variables['u'][:])
    v_vel = np.squeeze(fid.variables['v'][:])
    # Interpolate
    u_vel_2 = (u_vel[:,1:-1,:-1]+u_vel[:,1:-1,1:])/2
    v_vel_2 = (v_vel[:,:-1,1:-1]+v_vel[:,1:,1:-1])/2
    #u_vel_2 = (u_vel[:,:,:-1]+u_vel[:,:,1:])/2
    #v_vel_2 = (v_vel[:,:-1,:]+v_vel[:,1:,:])/2 
    print u_vel_2.shape
    print v_vel_2.shape
    # projection
    u4 = u_vel_2*(np.cos(angs)) + v_vel_2*(np.cos(np.pi/2 + angs))
    v4 = u_vel_2*(np.sin(angs)) + v_vel_2*(np.sin(np.pi/2 + angs))

    mag = np.sqrt((u_vel_2)**2+(v_vel_2)**2)
    store[:,1:-1,1:-1] = mag
    return store

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

grd = pyroms.grid.get_ROMS_grid('CCS')
#glat = grd.hgrid.lat_rho[1:-1,1:-1]
#glon = grd.hgrid.lon_rho[1:-1,1:-1]
angs = grd.hgrid.angle_rho[1:-1,1:-1]
x,y  = angs.shape 
print angs.shape
istart = int(3)
iend   = int(302)
jstart = int(486)
jend   = int(326)

salt = get_salt()
cs = get_vel()
print cs.shape
#cs_transect, z, lon, lat = pyroms.tools.transect(cs, istart, iend, jstart, jend, grd,vert=False, Cpos='rho')

#x = np.tile(np.array(range(transect[:].shape[1])),(transect[:].shape[0],1))
#v = np.linspace(-1*col_val,col_val, 50, endpoint=True)

#fig = plt.figure()
#ax = fig.add_subplot(111,axisbg=[0.5,0.5,0.5])
#v = np.linspace(-1*col_val,col_val, 50, endpoint=True)
#cs = ax.contourf(x,z,transect,v,extend = 'both', cmap='bwr')
plt.colorbar(cs,ticks=[-.001, 0, 0.001])
plt.show()

