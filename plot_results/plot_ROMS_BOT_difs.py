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

    # BOTTOM
    nc1 = '/glade/p/work/edrenkar/MODELS/CCS/ANALYSES/CCS-LD.HCo02Y/5yr_his_BOT_TEMP_clim.nc'
    nc2 = '/glade/p/work/edrenkar/MODELS/CCS/ANALYSES/CCS-LD.FCo017/5yr_fut_BOT_TEMP_clim.nc'

    fid1 = nc.Dataset(nc1)
    sst1 = fid1.variables['temp'][i,:].squeeze()
    
    fid2 = nc.Dataset(nc2)
    sst2 = fid2.variables['temp'][i,:].squeeze()
    sst = sst2-sst1 
    print np.min(sst),np.max(sst)
    return sst,mon

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ROMS Grid information
grdfile = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid = nc.Dataset(grdfile)
mask_rho = fid.variables['mask_rho'][:]
rlat = fid.variables['lat_rho'][:]
rlon = fid.variables['lon_rho'][:]
plat = fid.variables['lat_psi'][:]
plon = fid.variables['lon_psi'][:]
dep = fid.variables['h'][:]

### OFFSETS
joffset = 0
ioffset = 0

m_offset = 0.05
mask_val = 0
map_order = 30

# INITIAL FIGURE
#fig, ax = plt.subplots(figsize=(4,6))
fig, ax = plt.subplots(figsize=(4,8))
#fig, ax = plt.subplots(figsize=(4,4))
#fig, ax = plt.subplots(figsize=(4,6))
fig.tight_layout()
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

tx =plt.text(-113,31.5,'', fontsize=20,zorder=map_order+5)

# TOP PLOT
#m.drawmeridians([-128,-124], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
#m.drawparallels([42,46,50], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)

# TOP MIDDLE PLOT
m.drawmeridians([-128,-124], labels=[0,0,0,0], linewidth=2,fmt='%d', fontsize=18,zorder=map_order+5)
m.drawparallels([36,40], labels=[0,0,0,0], linewidth=2,fmt='%d', fontsize=18,zorder=map_order+5)

# BOTTOM MIDDLE PLOT
#m.drawmeridians([-120,-116], labels=[0,0,0,0], linewidth=2.5,fmt='%d', fontsize=18,zorder=map_order+5)
#m.drawparallels([32,36], labels=[0,0,0,0], linewidth=2.5,fmt='%d', fontsize=18,zorder=map_order+5)

# BOTTOM PLOT
#m.drawmeridians([-116,-112], labels=[0,0,0,0], linewidth=1,fmt='%d', fontsize=18,zorder=map_order+5)
#m.drawparallels([24,28,32], labels=[0,0,0,0], linewidth=1,fmt='%d', fontsize=18,zorder=map_order+5)

sst,mon = get_sst(0)
im1 = m.pcolormesh(plon,plat,sst[1:-1,1:-1],vmin=0,vmax=4,cmap='nipy_spectral',zorder=map_order)
im2 = m.contour(rlon,rlat,dep,[100,200],linewidths=1,linestyles='dotted',colors='white',zorder=map_order+1)
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
    im1   = m.pcolormesh(plon,plat,sst[1:-1,1:-1],vmin=0,vmax=4,cmap='nipy_spectral',zorder=map_order)
    im2   = m.contour(rlon,rlat,dep,[100,200],linewidths=.5,linestyles='dotted',colors='white',zorder=map_order+1)
    polygon_patch(m,ax)
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    tx_str = mon 
    tx.set_text(tx_str)
    # ADD Colorbar on first iteration 
    #if i == 0:
    #   cbar = m.colorbar(im1, location='bottom',size="5%", pad="3%",ticks=[])

# TOP PLOT
#ax.set_ylim(40,np.max(plat)+m_offset-.4)
#ax.set_xlim(-128.8,-123.8)

# TOP MIDDLE PLOT
ax.set_ylim(35.5,41)
ax.set_xlim(-125,-121.5)

# BOTTOM MIDDLE PLOT
#ax.set_ylim(31.5,35)
#ax.set_xlim(-121,-116.5)

# BOTTOM PLOT
#ax.set_ylim(24.4,32)
#ax.set_xlim(-117,-112+m_offset)

ani = animation.FuncAnimation(fig, updatefig,frames=12, blit=False)
# ani.save('CCS_ROMS_SST_DIFS.gif', writer = 'ImageMagickWriter',fps=1)
# ani.save('CCS_ROMS_SST_DIFS.gif', writer = 'imagemagick',fps=1)
ani.save('CCS_ROMS_BOT_DIFS.gif', writer = 'imagemagick',fps=1)
#ani.save('CCS_ROMS_VEL_CLIM.gif', writer = 'imagemagick',fps=1)
#plt.show()
