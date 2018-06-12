import pyroms
import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection

def polygon_patch(mapid,axs):
    mapid.drawcoastlines(linewidth=0,zorder=map_order+6)
    mapid.drawmapboundary(fill_color=[.9,.97,1])
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor=(1,1,1), closed=False)
    axs.add_collection(lc)

def outline_mask(mapid,mask_img,val,x0,y0,x1,y1):
    mapimg = (mask_img == val)
    ver_seg = np.where(mapimg[:,1:] != mapimg[:,:-1])
    hor_seg = np.where(mapimg[1:,:] != mapimg[:-1,:])

    l = []
    v = []
    # horizonal segments
    for p in zip(*hor_seg):
        if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
           v.append((glon[p[0]+1,p[1]],glat[p[0]+1,p[1]]))
        else :
           l.append((glon[p[0]+1,p[1]],glat[p[0]+1,p[1]]))

        if p[1] == mask_img.shape[1] - 1 :
           if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
               v.append((glon[p[0]+1,p[1]],glat[p[0]+1,p[1]]))
           else:
               l.append((glon[p[0]+1,p[1]],glat[p[0]+1,p[1]]))
        else :
           if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
              v.append((glon[p[0]+1,p[1]+1],glat[p[0]+1,p[1]+1]))
           else:
              l.append((glon[p[0]+1,p[1]+1],glat[p[0]+1,p[1]+1]))

        l.append((np.nan,np.nan))
        v.append((np.nan,np.nan))
    #vertical segments
    for p in zip(*ver_seg):
        if p[1] == mask_img.shape[1]-1:
           if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
              v.append((glon[p[0],p[1]],glat[p[0],p[1]]))
              v.append((glon[p[0]+1,p[1]],glat[p[0]+1,p[1]]))
           else:
              l.append((glon[p[0],p[1]],glat[p[0],p[1]]))
              l.append((glon[p[0]+1,p[1]],glat[p[0]+1,p[1]]))
        elif p[0] == mask_img.shape[0]-1:
             if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
              v.append((glon[p[0],p[1]],glat[p[0],p[1]]))
              v.append((glon[p[0]+1,p[1]],glat[p[0]+1,p[1]]))
             else:
              l.append((glon[p[0],p[1]],glat[p[0],p[1]]))
              l.append((glon[p[0],p[1]+1],glat[p[0],p[1]+1]))
        else:
           if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
              v.append((glon[p[0],p[1]+1],glat[p[0],p[1]+1]))
              v.append((glon[p[0]+1,p[1]+1],glat[p[0]+1,p[1]+1]))
           else:
              l.append((glon[p[0],p[1]+1],glat[p[0],p[1]+1]))
              l.append((glon[p[0]+1,p[1]+1],glat[p[0]+1,p[1]+1]))

        l.append((np.nan, np.nan))
        v.append((np.nan, np.nan))
    segments = np.array(l)
    vip_segments = np.array(v)
    mapid.plot(segments[:,0], segments[:,1], latlon=True, color=(0,0,0), linewidth=.75,zorder=map_order+2)
    mapid.plot(vip_segments[:,0], vip_segments[:,1], latlon=True, color=(0,0,0), linewidth=.75,zorder=map_order+3)

def fill_CA_Gulf(field):
    field[30:33,65:]=-1
    field[28:30,66:]=-1
    field[24:28,67:]=-1
    return field    

def get_sst(i):
    MONS = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    mon = MONS[i] 
    #ncfile = '/Users/elizabethdrenkard/Desktop/CESM_TEMP/016_TEMP_' + str(i+1).zfill(2) + '_delta.nc'
    ncfile = '/Users/elizabethdrenkard/Desktop/ECCWO_FILES/SST_DELTA_017.nc' 
    fid = nc.Dataset(ncfile)
    sst = np.ma.masked_array(fid.variables['temp'][i,0,:].squeeze())    
    sst[sst >100] = -1
    print sst.shape 
    sst = fill_CA_Gulf(sst)
    #lat = fid.variables['lat'][:] 
    #lon = fid.variables['lon'][:]
    #lon[lon>180]=lon[lon>180]-360   
    #return sst,lat,lon, mon
    return sst, mon
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# CCS grid shape ONLY
GRD = pyroms.grid.get_ROMS_grid('CCS')
mask = GRD.hgrid.mask_rho
glat = GRD.hgrid.lat_rho
glon = GRD.hgrid.lon_rho

### OFFSETS
joffset = 0
ioffset = 0

m_offset = 0.05
mask_val = 0
map_order = 30
vip_eta = [0,872]
vip_xi  = [240,378]

# INITIAL FIGURE
fig, ax = plt.subplots(figsize=(8,8))
m = Basemap(llcrnrlat=np.min(glat)-m_offset,urcrnrlat = np.max(glat)+m_offset,llcrnrlon=np.min(glon)-m_offset,urcrnrlon=np.max(glon)+m_offset, resolution='i', ax=ax)

#P = m.pcolormesh(glon,glat,mask,vmin=.5,vmax=.75,edgecolors='face',cmap='Blues',zorder=map_order)
#P.cmap.set_under('white')
#P.cmap.set_over([.9,.97,1])

grd_fil = '/Users/elizabethdrenkard/Desktop/ECCWO_FILES/LENS_grid.nc'
lsm = nc.Dataset(grd_fil).variables['lsm'][:]
lats = nc.Dataset(grd_fil).variables['lat'][:]
lons = nc.Dataset(grd_fil).variables['lon'][:]
lons[lons>180]=lons[lons>180]-360
outline_mask(m,mask,mask_val,glon[0,0],glat[0,0],glon[-1,-1],glat[-1,-1])

#DOMAIN OUTLINE
for j in range(glat.shape[0]-2):
    m.plot((glon[j,0],glon[j+1,0]),(glat[j,0],glat[j+1,0]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((glon[j,-1],glon[j+1,-1]),(glat[j,-1],glat[j+1,-1]),linewidth=2,color='k',zorder=map_order+1)
for ii in range(glat.shape[1]-2):
    m.plot((glon[0,ii],glon[0,ii+1]),(glat[0,ii],glat[0,ii+1]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((glon[-1,ii],glon[-1,ii+1]),(glat[-1,ii],glat[-1,ii+1]),linewidth=2,color='k',zorder=map_order+1)

#plt.title('SST (oC)')
tx =plt.text(-118,46.5,'', fontsize=20,zorder=map_order+5)
n=0

m.drawmeridians([-142,-111], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
m.drawparallels([18,50], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)

sst,mon = get_sst(0)
im1 = m.pcolormesh(lons,lats,sst[:],vmin=2,vmax=5,cmap='nipy_spectral',zorder=map_order)
im1.cmap.set_under('w')
polygon_patch(m,ax)
# ANIMATION
def updatefig(i):
    global im1,im2,tx
    print i
    # REMOVE images after first step
    if i > 0:
       im1.remove()
    sst, mon = get_sst(i)
    im1 = m.pcolormesh(lons,lats,sst[:],vmin=2,vmax=5,cmap='nipy_spectral',zorder=map_order)
    im1.cmap.set_under('w')
    polygon_patch(m,ax)
    tx_str = mon 
    tx.set_text(tx_str)
   
    # ADD Colorbar on first iteration 
    if i == 0:
       cbar = m.colorbar(im1, location='bottom',size="5%", pad="3%",ticks=[])

ani = animation.FuncAnimation(fig, updatefig,frames=12, blit=False)
ani.save('CESM_SST_DELTA_CLIM.gif', writer = 'imagemagick',fps=1)
plt.show()
