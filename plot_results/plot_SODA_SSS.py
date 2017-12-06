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

def fill_CA_Gulf(field):
    A = 180*np.ones(4,dtype=np.int)
    B = 178*np.ones(1,dtype=np.int)
    BB = 176*np.ones(2,dtype=np.int)
    C = 174*np.ones(4,dtype=np.int)
    D = 171*np.ones(8,dtype=np.int)
    E = 164*np.ones(11,dtype=np.int)
    F = 159*np.ones(12,dtype=np.int)
    nc = np.append(A,np.append(B,np.append(BB,np.append(C,np.append(D,np.append(E,F))))))
    for nr in range(40,80+1):
        field[nr,nc[nr-40]:]=ma.masked
        
    return field    

def get_sss(ncfil):
    fid = nc.Dataset(ncfil)
    sst = fid.variables['salt'][0,0,:].squeeze()    
    print sst.shape
    sst = fill_CA_Gulf(sst)
    lat = fid.variables['yt_ocean'][:] 
    lon = fid.variables['xt_ocean'][:]   
    return sst,lat,lon

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# DATA LOCATION 
ncfil = '/glade/p/work/edrenkar/external_data/SODA/soda_annual_avg_salt.nc'
sss, lat, lon = get_sss(ncfil)

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

clon,clat = np.meshgrid(lon,lat)

# INITIAL FIGURE
fig, ax = plt.subplots(figsize=(8,8))
m = Basemap(llcrnrlat=np.min(glat)-m_offset,urcrnrlat = np.max(glat)+m_offset,llcrnrlon=np.min(glon)-m_offset,urcrnrlon=np.max(glon)+m_offset, resolution='i', ax=ax)

P = m.contourf(clon,clat,sss,np.linspace(30,36+1,50),edgecolors='face',cmap='plasma',zorder=map_order)
P.cmap.set_under('white')
P.cmap.set_over([.9,.97,1])
plt.colorbar(P)
##DOMAIN OUTLINE
#for j in range(glat.shape[0]-2):
#    m.plot((glon[j,0],glon[j+1,0]),(glat[j,0],glat[j+1,0]),linewidth=2,color='k',zorder=map_order+1)
#    m.plot((glon[j,-1],glon[j+1,-1]),(glat[j,-1],glat[j+1,-1]),linewidth=2,color='k',zorder=map_order+1)
#for ii in range(glat.shape[1]-2):
#    m.plot((glon[0,ii],glon[0,ii+1]),(glat[0,ii],glat[0,ii+1]),linewidth=2,color='k',zorder=map_order+1)
#    m.plot((glon[-1,ii],glon[-1,ii+1]),(glat[-1,ii],glat[-1,ii+1]),linewidth=2,color='k',zorder=map_order+1)

polygon_patch(m,ax)

m.drawmeridians([-142,-111], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
m.drawparallels([18,50], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)

plt.show()
