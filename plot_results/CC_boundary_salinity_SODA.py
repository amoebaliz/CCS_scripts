import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import colors
import pyroms
from matplotlib.collections import PolyCollection
from mpl_toolkits.basemap import Basemap

def polygon_patch(mapid,axs):
    mapid.drawcoastlines(linewidth=0)
    mapid.drawmapboundary(fill_color=[.9,.97,1])
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor=(0,0,0), closed=False)
    axs.add_collection(lc)

# ~~~~~~~~~~~~~

# MODEL DATA FILES 
mean_salt_fil = '/glade/p/work/edrenkar/external_data/SODA/soda_annual_avg_salt.nc'

# SALT DATA
fid = nc.Dataset(mean_salt_fil)
salt = fid.variables['salt'][:].squeeze()

#SODA GRID DATA
lon1d = fid.variables['xt_ocean'][:].squeeze()
lat1d = fid.variables['yt_ocean'][:].squeeze()
lon,lat = np.meshgrid(lon1d,lat1d)
z = fid.variables['st_ocean'][:].squeeze()

# Limiting SODA TO UPPER 1000m
IZ = np.where(z<1000)[0]

# FIGURE DETAILS
cmap_file = 'auad_salinity.cmap'
cmap_a = np.loadtxt(cmap_file)/256.
cmap = colors.ListedColormap(cmap_a)
m_offset = 0
mask_val = 0
map_order = 30

# Find and Store salinity minimum
smin = np.amin(salt[IZ,:,:],axis=(2,0))

# Initialize variable lat_x_lon variable to store depths
depth_store = np.zeros(salt.shape[1:3])

# MAP FIGURE
fig = plt.figure(figsize=(10,10))
fig.subplots_adjust(left=.1, right=.9, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False)

# BASEMAP OBJECT
m = Basemap(llcrnrlat=20-m_offset,urcrnrlat = 50+m_offset, llcrnrlon=-150-m_offset, urcrnrlon=-110+m_offset, resolution='f', ax=ax)


# Identify Regions that are within 0.5 psu of minimum
for ny in range(salt.shape[1]):
    depth_slab = ma.masked_array(np.zeros((len(IZ),salt.shape[2])))
    tmp_salt = salt[IZ,ny,:].squeeze()

    # SODA EVALUATION
    #Initial min salt eval
    smin = np.min(tmp_salt)
    if ((ny>100) & (ny<=140)):
       smin = np.min(tmp_salt[:,:(tmp_salt[0,:].count()-6)]) 
       smin = np.max((smin,32))
    elif ny > 140:
       smin = np.max((smin,32))

    print ny,smin
    #if ((ny > 100)) & ((ny<=137)):
    #   smin = np.min(tmp_salt[:,:(tmp_salt[0,:].count()-6)]) 
    #elif ((ny > 137)&(ny<182)):
    #   smin = np.min(tmp_salt[:,:(tmp_salt[0,:].count()-8)])
    #   smin = 32
    #print ny, smin
    for nc in range(salt.shape[2]):
        depth_col = ma.masked_array(np.zeros(len(IZ)))
        depth_col[tmp_salt[:,nc].squeeze() <= (smin+0.5)] = z[tmp_salt[:,nc].squeeze() <= (smin+0.5)]
# SURFACE CRITERIA: SODA
        if ((depth_col[0] == 0)):# or (depth_col[:])):   
           depth_col[:] = 0
        depth_slab[:,nc] = depth_col
# MASK CRITERIAL: SODA
    depth_slab[ma.getmask(tmp_salt)]=0   

    depth_store[ny,:] = np.amax(depth_slab,axis=0)
depth_store[depth_store==0]=np.nan

#SODA
P = m.contourf(lon,lat,depth_store,10*np.arange(2,28,2),edgecolors='face',cmap=cmap,zorder=map_order)
C = m.contour(lon,lat,depth_store,[100,200], colors='k',zorder=map_order)

cbar_ax = fig.add_axes([0.91,0.15,0.02,0.7])
plt.colorbar(P,cax=cbar_ax)

polygon_patch(m,ax)

m.drawmeridians([-150,-110], labels=[0,0,1,0], fmt='%d', fontsize=18,zorder=map_order+5)
m.drawparallels([20,50], labels=[1,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
#m.drawmeridians([-142,-111], labels=[0,0,1,0], fmt='%d', fontsize=18,zorder=map_order+5)
#m.drawparallels([18,50], labels=[1,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
plt.show()

