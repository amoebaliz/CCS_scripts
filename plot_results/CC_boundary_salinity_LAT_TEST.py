import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import colors
import pyroms
import ESMF
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

def outline_mask(mapid,mask_img,val,x0,y0,x1,y1):
    mapimg = (mask_img == val)
    ver_seg = np.where(mapimg[:,1:] != mapimg[:,:-1])
    hor_seg = np.where(mapimg[1:,:] != mapimg[:-1,:])

    l = []
    v = []
    # horizonal segments
    for p in zip(*hor_seg):
        if (ccs_eta[0] < p[0] < ccs_eta[1] and ccs_xi[0] < p[1] < ccs_xi[1]):
           v.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))
        else :
           l.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))

        if p[1] == mask_img.shape[1] - 1 :
           if (ccs_eta[0] < p[0] < ccs_eta[1] and ccs_xi[0] < p[1] < ccs_xi[1]):
               v.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))
           else:
               l.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))
        else :
           if (ccs_eta[0] < p[0] < ccs_eta[1] and ccs_xi[0] < p[1] < ccs_xi[1]):
              v.append((lon[p[0]+1,p[1]+1],lat[p[0]+1,p[1]+1]))
           else:
              l.append((lon[p[0]+1,p[1]+1],lat[p[0]+1,p[1]+1]))

        l.append((np.nan,np.nan))
        v.append((np.nan,np.nan))
    #vertical segments
    for p in zip(*ver_seg):
        if p[1] == mask_img.shape[1]-1:
           if (ccs_eta[0] < p[0] < ccs_eta[1] and ccs_xi[0] < p[1] < ccs_xi[1]):
              v.append((lon[p[0],p[1]],lat[p[0],p[1]]))
              v.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))
           else:
              l.append((lon[p[0],p[1]],lat[p[0],p[1]]))
              l.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))
        elif p[0] == mask_img.shape[0]-1:
             if (ccs_eta[0] < p[0] < ccs_eta[1] and ccs_xi[0] < p[1] < ccs_xi[1]):
              v.append((lon[p[0],p[1]],lat[p[0],p[1]]))
              v.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))
             else:
              l.append((lon[p[0],p[1]],lat[p[0],p[1]]))
              l.append((lon[p[0],p[1]+1],lat[p[0],p[1]+1]))
        else:
           if (ccs_eta[0] < p[0] < ccs_eta[1] and ccs_xi[0] < p[1] < ccs_xi[1]):
              v.append((lon[p[0],p[1]+1],lat[p[0],p[1]+1]))
              v.append((lon[p[0]+1,p[1]+1],lat[p[0]+1,p[1]+1]))
           else:
              l.append((lon[p[0],p[1]+1],lat[p[0],p[1]+1]))
              l.append((lon[p[0]+1,p[1]+1],lat[p[0]+1,p[1]+1]))

        l.append((np.nan, np.nan))
        v.append((np.nan, np.nan))
    segments = np.array(l)
    ccs_segments = np.array(v)
    mapid.plot(segments[:,0], segments[:,1], latlon=True, color=(0,0,0), linewidth=.75,zorder=map_order+2)
    mapid.plot(ccs_segments[:,0], ccs_segments[:,1], latlon=True, color=(0,0,0), linewidth=.75,zorder=map_order+3)

# ~~~~~~~~~~~~~

# MODEL DATA FILES 
mean_salt_fil = '/Volumes/Abalone/CCS/his2/CCS-LD.HCo02Y_10yr_clim_salt.nc'

# SALT DATA
fid = nc.Dataset(mean_salt_fil)
salt = np.mean(fid.variables['salt'][:],axis=0).squeeze()

## CCS ROMS GRID INFORMATION
grd = pyroms.grid.get_ROMS_grid('CCS')
mask = grd.hgrid.mask_rho
lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho
h = grd.vgrid.h[:]
hc = grd.vgrid.hc
N = grd.vgrid.N
s_rho = grd.vgrid.s_rho[:]
Cs_r = grd.vgrid.Cs_r[:]
Vtrans = grd.vgrid.Vtrans
zeta = np.mean(fid.variables['zeta'][:],axis=0).squeeze()
z = pyroms.vgrid.z_r(h,hc,N,s_rho,Cs_r,zeta,Vtrans)[:]

# FIGURE DETAILS
cmap_file = 'auad_salinity.cmap'
cmap_a = np.loadtxt(cmap_file)/256.
cmap = colors.ListedColormap(cmap_a)
m_offset = 0
mask_val = 0
map_order = 30
ccs_eta = [0,872]
ccs_xi  = [240,378]

# Initialize variable lat_x_lon variable to store depths
depth_store = np.ma.ones(salt.shape[1:3])*10
# MAP FIGURE
fig = plt.figure(figsize=(10,10))
fig.subplots_adjust(left=.1, right=.9, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False)

# BASEMAP OBJECT
m = Basemap(llcrnrlat=20-m_offset,urcrnrlat = 50+m_offset, llcrnrlon=-150-m_offset, urcrnrlon=-110+m_offset, resolution='f', ax=ax)

lat_inc = 0.25
LATS = np.arange(np.ceil(np.max(lat)), np.floor(np.min(lat))-1,-1*lat_inc)
for nlat in LATS:
    lat_mask = np.zeros(lat.shape)
    lat_mask[(lat<nlat-lat_inc/2)]=1
    lat_mask[(lat>nlat+lat_inc/2)]=1

    tmp_salt = salt.copy()
    tmp_salt[:,lat_mask.astype(np.bool)]=np.ma.masked 

    # salt min at given latitude
    smin = np.ma.min(tmp_salt)

    tmp_z = np.ma.masked_array(z.copy())

    # MASK ANYTHING LAND AND NOT IN LAT RANGE
    tmp_z[tmp_salt.mask] = np.ma.masked

    # MASK ANYTHIN OVER .5+smin
    tmp_z[tmp_salt > (smin + 0.5)]  = np.ma.masked  
    
    # MASK ANYTHING WHERE MIN IS NOT AT SURFACE 
    #tmp_z[:,tmp_z[-1,:].mask] = np.ma.masked
    print nlat, smin

    #if (tmp_z[-1,:].all() is np.ma.masked):
    #   break

    in_lat_mask = -1*lat_mask+1
 
    depth_store[in_lat_mask.astype(np.bool)] = np.ma.min(tmp_z,axis=0)[in_lat_mask.astype(np.bool)]
    depth_store[depth_store>5]=np.ma.masked
    #plt.pcolor(lon,lat,np.ma.max(tmp_z,axis=0).squeeze(),vmin=-1,vmax=1)
#plt.pcolor(lon,lat,depth_store,vmin=-500,vmax=0)
#plt.colorbar()
#plt.show()


P = m.contourf(lon,lat,depth_store*-1,10*np.arange(2,28,2),edgecolors='face',cmap=cmap,zorder=map_order)
C = m.contour(lon,lat,depth_store*-1,[100,200], colors='k',zorder=map_order)

cbar_ax = fig.add_axes([0.91,0.15,0.02,0.7])
plt.colorbar(P,cax=cbar_ax)

# ROMS DOMAIN OUTLINE
outline_mask(m,mask,mask_val,lon[0,0],lat[0,0],lon[-1,-1],lat[-1,-1])

for j in range(lat.shape[0]-2):
    m.plot((lon[j,0],lon[j+1,0]),(lat[j,0],lat[j+1,0]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((lon[j,-1],lon[j+1,-1]),(lat[j,-1],lat[j+1,-1]),linewidth=2,color='k',zorder=map_order+1)
for ii in range(lat.shape[1]-2):
    m.plot((lon[0,ii],lon[0,ii+1]),(lat[0,ii],lat[0,ii+1]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((lon[-1,ii],lon[-1,ii+1]),(lat[-1,ii],lat[-1,ii+1]),linewidth=2,color='k',zorder=map_order+1)

# SUBDOMAIN OUTLINE
for j in range(200,700):
    m.plot((lon[j,0],lon[j+1,0]),(lat[j,0],lat[j+1,0]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((lon[j,-1],lon[j+1,-1]),(lat[j,-1],lat[j+1,-1]),linewidth=2,color='k',zorder=map_order+1)
for ii in range(79,378):
    m.plot((lon[0,ii],lon[0,ii+1]),(lat[0,ii],lat[0,ii+1]),linewidth=1,color='k',zorder=map_order+2)
    m.plot((lon[-1,ii],lon[-1,ii+1]),(lat[-1,ii],lat[-1,ii+1]),linewidth=1,color='k',zorder=map_order+2)

polygon_patch(m,ax)

m.drawmeridians([-150,-110], labels=[0,0,1,0], fmt='%d', fontsize=18,zorder=map_order+5)
m.drawparallels([20,50], labels=[1,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
#plt.savefig('10yr_his_salt_min.png')
plt.show()

