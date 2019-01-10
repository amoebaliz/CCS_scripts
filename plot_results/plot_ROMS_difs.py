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

def get_sst(i):
    MONS = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    mon = MONS[i]
 
    # Surface
    #nc1 = '/Volumes/Abalone/CCS/his/clim/SST_10y_clim.nc'
    #nc2 = '/Volumes/Abalone/CCS/016/clim/SST_10y_clim.nc'
    #nc1 = '/Users/elizabethdrenkard/Desktop/ECCWO_FILES/SST_5yr_his_clim.nc'    
    #nc2 = '/Users/elizabethdrenkard/Desktop/ECCWO_FILES/SST_5yr_fut_clim.nc'  

    nc1 = '/Users/elizabethdrenkard/Desktop/SST_10yr_his_clim.nc'    
    nc2 = '/Users/elizabethdrenkard/Desktop/SST_10yr_fut017_clim.nc'  

    # 50-meter depth
    # nc1 = '/Volumes/Abalone/CCS/his/clim/T50_10y_clim.nc'
    # nc2 = '/Volumes/Abalone/CCS/016/clim/T50_10y_clim.nc'

    # BOTTOM
    #nc1 = '/Volumes/Abalone/CCS/his/clim/BOT_10y_clim.nc'
    #nc2 = '/Volumes/Abalone/CCS/016/clim/BOT_10y_clim.nc'


    fid1 = nc.Dataset(nc1)
    sst1 = fid1.variables['temp'][i,:].squeeze()
    
    fid2 = nc.Dataset(nc2)
    sst2 = fid2.variables['temp'][i,:].squeeze()
    sst = sst2-sst1 
    #print np.min(sst),np.max(sst)
    return sst,mon


def get_vel(i):
    MONS = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    mon = MONS[i]
    nc1 = '/Users/elizabethdrenkard/Desktop/CCS_ROMS_CLIM_FLS/CCS_5y_his_UV_clim.nc'
    nc2 = '/Users/elizabethdrenkard/Desktop/CCS_ROMS_fut/Fut_2yr_'+str(i+1).zfill(2)+'.nc'

    fid1 = nc.Dataset(nc1)
    fid2 = nc.Dataset(nc2)

    u_vel1 = np.squeeze(fid1.variables['u'][i,:])
    v_vel1 = np.squeeze(fid1.variables['v'][i,:])

    u_vel2 = np.squeeze(fid2.variables['u'][0,:])
    v_vel2 = np.squeeze(fid2.variables['v'][0,:])

    u_vel = u_vel2-u_vel1
    v_vel = v_vel2-v_vel2

    # Interpolate
    u_vel_2 = (u_vel[1:-1,:-1]+u_vel[1:-1,1:])/2
    v_vel_2 = (v_vel[:-1,1:-1]+v_vel[1:,1:-1])/2
    # projection
    u4 = u_vel_2*(np.cos(angs)) + v_vel_2*(np.cos(np.pi/2 + angs))
    v4 = u_vel_2*(np.sin(angs)) + v_vel_2*(np.sin(np.pi/2 + angs))

    mag = np.sqrt((u_vel_2)**2+(v_vel_2)**2)
    return u4, v4, mag, mon

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ROMS Grid information
grdfile = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid = nc.Dataset(grdfile)
mask_rho = fid.variables['mask_rho'][:]
rlat = fid.variables['lat_rho'][:]
rlon = fid.variables['lon_rho'][:]
plat = fid.variables['lat_psi'][:]
plon = fid.variables['lon_psi'][:]

### OFFSETS
joffset = 0
ioffset = 0

m_offset = 0.05
mask_val = 0
map_order = 30

# INITIAL FIGURE
fig, ax = plt.subplots(figsize=(8,8))
m = Basemap(llcrnrlat=np.min(plat)-m_offset,urcrnrlat = np.max(plat)+m_offset,llcrnrlon=np.min(plon)-m_offset,urcrnrlon=np.max(plon)+m_offset, resolution='i', ax=ax)

P = m.pcolormesh(plon,plat,mask_rho[1:-1,1:-1],vmin=.5,vmax=.75,edgecolors='face',cmap='Blues',zorder=map_order)
P.cmap.set_under('white')
P.cmap.set_over([.9,.97,1])

outline_mask(m,mask_rho[1:-1,1:-1],mask_val,plon[0,0],plat[0,0],plon[-1,-1],plat[-1,-1])

#DOMAIN OUTLINE
for j in range(plat.shape[0]-1):
    m.plot((plon[j,0],plon[j+1,0]),(plat[j,0],plat[j+1,0]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((plon[j,-1],plon[j+1,-1]),(plat[j,-1],plat[j+1,-1]),linewidth=2,color='k',zorder=map_order+1)
for ii in range(plat.shape[1]-1):
    m.plot((plon[0,ii],plon[0,ii+1]),(plat[0,ii],plat[0,ii+1]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((plon[-1,ii],plon[-1,ii+1]),(plat[-1,ii],plat[-1,ii+1]),linewidth=2,color='k',zorder=map_order+1)

polygon_patch(m,ax)
afreq = 25
#plt.title('SST (oC)')
tx =plt.text(-118,46.5,'', fontsize=20,zorder=map_order+5)
n=0

m.drawmeridians([-142,-111], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
m.drawparallels([18,50], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)

#u,v,mag,mon = get_vel(0)
sst,mon = get_sst(0)
#im1 = m.pcolor(plon,plat,sst[1:-1,1:-1],vmin=-5,vmax=5,cmap='bwr',zorder=map_order)
im1 = m.pcolor(plon,plat,sst[1:-1,1:-1],vmin=2,vmax=5,cmap='nipy_spectral',zorder=map_order)
#u,v,mag,mon = get_vel(0)
#im1 = m.pcolor(plon,plat,mag,vmin=0,vmax=.2,cmap='OrRd',zorder=map_order)
#im2 = m.quiver(rlon[::afreq],rlat[::afreq],u[::afreq,::afreq],v[::afreq,::afreq], scale=5,zorder=map_order+2)

#plt.show()
# ANIMATION
def updatefig(i):
    global im1,im2,tx
    print i
    # REMOVE images after first step
    #im1.remove()
    # im2.remove()
    #im2.remove()
    if i > 0:
       im1.remove()
    sst,mon = get_sst(i)
    print np.mean(sst), np.std(sst)
    #u,v,mag,mon = get_vel(i)
    #im1   = m.pcolormesh(plon,plat,sst[1:-1,1:-1],vmin=-5,vmax=5,cmap='bwr',zorder=map_order)
    im1   = m.pcolormesh(plon,plat,sst[1:-1,1:-1],vmin=2,vmax=5,cmap='nipy_spectral',zorder=map_order)
    polygon_patch(m,ax)
    tx_str = mon 
    tx.set_text(tx_str)
    # ADD Colorbar on first iteration 
    if i == 0:
       cbar = m.colorbar(im1, location='bottom',size="5%", pad="3%",ticks=[])

ani = animation.FuncAnimation(fig, updatefig,frames=12, blit=False)
ani.save('CCS_ROMS_SST_DIFS.gif', writer = 'imagemagick',fps=1)
plt.show()
