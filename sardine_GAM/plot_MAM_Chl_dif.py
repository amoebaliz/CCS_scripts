import numpy as np
import numpy.ma as ma
import ESMF
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
        v.append((plon[p[0]+1,p[1]],plat[p[0]+1,p[1]]))
        v.append((plon[p[0]+1,p[1]+1],plat[p[0]+1,p[1]+1]))

        l.append((np.nan,np.nan))
        v.append((np.nan,np.nan))
    #vertical segments
    for p in zip(*ver_seg):
        l.append((plon[p[0],p[1]+1],plat[p[0],p[1]+1]))
        l.append((plon[p[0]+1,p[1]+1],plat[p[0]+1,p[1]+1]))

        l.append((np.nan, np.nan))
        v.append((np.nan, np.nan))

    l_segments = np.array(l)
    v_segments = np.array(v)
    mapid.plot(l_segments[:,0], l_segments[:,1], latlon=True, color=(0,0,0), linewidth=.75,zorder=map_order+2)
    mapid.plot(v_segments[:,0], v_segments[:,1], latlon=True, color=(0,0,0), linewidth=.75,zorder=map_order+3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ROMS Grid information
grdfile = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid = nc.Dataset(grdfile)
mask_rho = fid.variables['mask_rho'][:]
rlat = fid.variables['lat_rho'][:]
rlon = fid.variables['lon_rho'][:]
plat = fid.variables['lat_psi'][:]
plon = fid.variables['lon_psi'][:]

roms_mask = fid.variables['mask_rho'][:]
roms_mask = -1*roms_mask+1

ny = rlat.shape[0]
nx = rlon.shape[1]

# Chl INFO
#nc1 = 'SeaWiFS_CCS_his_Clim_MAM.nc'
#nc2 = 'SeaWiFS_CCS_fut_Clim_MAM.nc'

nc1 = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/CESM_LENS/CHL/CESMfut017_CCS_Clim_MAM.nc.full_rho'

fid1 = nc.Dataset(nc1)
chl = ma.array(fid1.variables['Chla'][:].squeeze())
chl[chl==0]=ma.masked

lon = fid1.variables['lon'][:]
lat = fid1.variables['lat'][:]

chl[chl==0] = ma.masked

#max_change = np.max(np.abs(chl_dif))
#his_range = np.max(chl1) - np.min(chl1)
#print np.max(chl1)
#print max_change
#print his_range

#print 100*max_change/his_range

# MAP SPECS
m_offset = 0.05
mask_val = 0
map_order = 30

### SST DIFERENCE FIGURE ######################################
fig, ax = plt.subplots(figsize=(8,8))
m = Basemap(llcrnrlat=np.min(plat)-m_offset,urcrnrlat = np.max(plat)+m_offset,llcrnrlon=np.min(plon)-m_offset,urcrnrlon=np.max(plon)+m_offset, resolution='i', ax=ax)

outline_mask(m,mask_rho[1:-1,1:-1],mask_val,plon[0,0],plat[0,0],plon[-1,-1],plat[-1,-1])

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

#im3 = m.pcolormesh(plon,plat,chl_dif,vmin=-.01,vmax=.01,cmap='PRGn',zorder=map_order)
#cbar = m.colorbar(im3, location='bottom',size="5%", pad="3%",ticks=[2,2.5,3,3.5,4])
#m.contour(rlon[1:-1,1:-1],rlat[1:-1,1:-1],chl_dif,0,linestyles='dashed',zorder=map_order+6)

print plon.shape, plat.shape, chl.shape

im3 = m.pcolormesh(plon,plat,chl[1:-1,1:-1],vmin=-.01,vmax=.01,cmap='PRGn',zorder=map_order)
cbar = m.colorbar(im3, location='bottom',size="5%", pad="3%",ticks=[2,2.5,3,3.5,4])
m.contour(rlon[1:-1,1:-1],rlat[1:-1,1:-1],chl[1:-1,1:-1],0,linestyles='dashed',zorder=map_order+6)

print np.min(chl)
plt.savefig('/Users/elizabethdrenkard/Desktop/Delta_MAM_Chl.png')
plt.show()
