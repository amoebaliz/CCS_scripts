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
#ncfil = '/glade/p/work/edrenkar/external_data/SODA/soda_annual_avg_salt.nc'
#ncfil = '/Users/elizabethdrenkard/Desktop/soda_annual_avg_salt.nc'
sss, lat, lon = get_sss(ncfil)

# ROMS Grid information
grdfile = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid = nc.Dataset(grdfile)
plat = fid.variables['lat_psi'][:]
plon = fid.variables['lon_psi'][:]


### OFFSETS
joffset = 0
ioffset = 0

m_offset = 0.05
map_order = 30

clon,clat = np.meshgrid(lon,lat)

# INITIAL FIGURE
fig, ax = plt.subplots(figsize=(8,8))
m = Basemap(llcrnrlat=np.min(plat)-m_offset,urcrnrlat = np.max(plat)+m_offset,llcrnrlon=np.min(plon)-m_offset,urcrnrlon=np.max(plon)+m_offset, resolution='i', ax=ax)

P = m.contourf(clon,clat,sss,np.linspace(32,35.5,50),cmap='inferno',zorder=map_order)
P.cmap.set_under('white')
P.cmap.set_over([.9,.97,1])
plt.colorbar(P)

#DOMAIN OUTLINE
for j in range(plat.shape[0]-1):
    m.plot((plon[j,0],plon[j+1,0]),(plat[j,0],plat[j+1,0]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((plon[j,-1],plon[j+1,-1]),(plat[j,-1],plat[j+1,-1]),linewidth=2,color='k',zorder=map_order+1)
for ii in range(plat.shape[1]-1):
    m.plot((plon[0,ii],plon[0,ii+1]),(plat[0,ii],plat[0,ii+1]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((plon[-1,ii],plon[-1,ii+1]),(plat[-1,ii],plat[-1,ii+1]),linewidth=2,color='k',zorder=map_order+1)

polygon_patch(m,ax)

m.drawmeridians([-142,-111], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
m.drawparallels([18,50], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)

plt.show()
