import numpy as np
import netCDF4 as nc
import pyroms
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection

def polygon_patch(mapid,axs):
    mapid.drawcoastlines(linewidth=0)
    mapid.drawmapboundary(fill_color=[.9,.97,1])
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor=(1,1,1), closed=False)
    axs.add_collection(lc)


grd = pyroms.grid.get_ROMS_grid('CCS')
lat = grd.hgrid.lat_rho[:]
lon = grd.hgrid.lon_rho[:]

fid = nc.Dataset('/Users/elizabethdrenkard/TOOLS/CCS_scripts/CESM_LENS/SSS_CESM_017_delta.nc')
salt = fid.variables['SSS'][:]
slat = fid.variables['lat'][:]
slon = fid.variables['lon'][:]

fid2 = nc.Dataset('/Users/elizabethdrenkard/external_data/CESM/drowned/drowned_LENS_017_SSS_delta.nc')
salt2 = fid2.variables['SSS'][:]
slat2 = fid2.variables['lat'][:]
slon2 = fid2.variables['lon'][:]
slon2[slon2>180]=slon2[slon2>180]-360

print slon2[0,:]

fig = plt.figure(figsize=(8,8))
fig.subplots_adjust(left=.1, right=.9, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False) #, xlim=(0, mask.shape[1]), ylim=(0, mask.shape[0]))

m_offset = 10
map_order = 3
m = Basemap(llcrnrlat=np.min(lat)-m_offset,urcrnrlat = np.max(lat)+m_offset,llcrnrlon=np.min(lon)-m_offset,urcrnrlon=np.max(lon)+m_offset, resolution='i', ax=ax)

P2 = m.pcolor(slon2,slat2,salt2[0,:].squeeze())

#DOMAIN OUTLINE
for j in range(lat.shape[0]-2):
    m.plot((lon[j,0],lon[j+1,0]),(lat[j,0],lat[j+1,0]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((lon[j,-1],lon[j+1,-1]),(lat[j,-1],lat[j+1,-1]),linewidth=2,color='k',zorder=map_order+1)
for ii in range(lat.shape[1]-2):
    m.plot((lon[0,ii],lon[0,ii+1]),(lat[0,ii],lat[0,ii+1]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((lon[-1,ii],lon[-1,ii+1]),(lat[-1,ii],lat[-1,ii+1]),linewidth=2,color='k',zorder=map_order+1)

polygon_patch(m,ax)
#plt.colorbar(P)


plt.show()
