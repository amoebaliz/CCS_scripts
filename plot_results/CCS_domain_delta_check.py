import pyroms
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt 
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# OBJECTIVE: determine whether deltas for Ocean and Atmosphere encompass CCS domain.

# CCS grid 
GRD = pyroms.grid.get_ROMS_grid('CCS')
mask = GRD.hgrid.mask_rho[:]
lat = GRD.hgrid.lat_rho[:]
lon = GRD.hgrid.lon_rho[:]

### OFFSETS
joffset = 0
ioffset = 0

m_offset =25
mask_val = 0
map_order = 30
ccs_eta = [0,870]
ccs_xi  = [0,376]


# INITIAL FIGURE
fig, ax = plt.subplots(figsize=(8,8))
m = Basemap(llcrnrlat=np.min(lat)-m_offset,urcrnrlat = np.max(lat)+m_offset,llcrnrlon=np.min(lon)-m_offset,urcrnrlon=np.max(lon)+m_offset, resolution='i', ax=ax)

# SWFSC FILES
ncfil = '/Users/liz.drenkard/Desktop/016_ocn_01.nc'
ncfil = '/Users/liz.drenkard/Desktop/LENS_016_msl_delta.nc'
ncfil = '/Users/liz.drenkard/Desktop/LENS_grid.nc'

fid = nc.Dataset(ncfil)
Slat = fid.variables['geolat_t'][:]
Slon = fid.variables['geolon_t'][:]
#Slon[Slon>180] = Slon[Slon>180]-360
Slon = Slon - 360
SSS  = fid.variables['ht'][:]

P = m.pcolormesh(Slon,Slat,SSS,zorder=map_order+1,cmap='spectral')
P.cmap.set_under('white')
outline_mask(m,mask,mask_val,lon[0,0],lat[0,0],lon[-1,-1],lat[-1,-1])

print np.min(Slon)
print np.max(Slon)



#DOMAIN OUTLINE
for j in range(lat.shape[0]-2):
    m.plot((lon[j,0],lon[j+1,0]),(lat[j,0],lat[j+1,0]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((lon[j,-1],lon[j+1,-1]),(lat[j,-1],lat[j+1,-1]),linewidth=2,color='k',zorder=map_order+1)
for ii in range(lat.shape[1]-2):
    m.plot((lon[0,ii],lon[0,ii+1]),(lat[0,ii],lat[0,ii+1]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((lon[-1,ii],lon[-1,ii+1]),(lat[-1,ii],lat[-1,ii+1]),linewidth=2,color='k',zorder=map_order+1)

polygon_patch(m,ax)
n=0

m.drawmeridians([-142,-111], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
m.drawparallels([18,50], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
plt.show()


