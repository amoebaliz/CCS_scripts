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
    mapid.plot(segments[:,0], segments[:,1], latlon=True, color=(0,0,0), linewidth=1,zorder=map_order+2)
    mapid.plot(vip_segments[:,0], vip_segments[:,1], latlon=True, color=(0,0,0), linewidth=1,zorder=map_order+3)


def get_sst(i):
    MONS = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    mon = MONS[i]

    # BOTTOM
    nc1 = '/Volumes/Abalone/CCS/his/clim/BOT_10y_clim.nc'
    nc2 = '/Volumes/Abalone/CCS/016/clim/BOT_10y_clim.nc'


    fid1 = nc.Dataset(nc1)
    sst1 = fid1.variables['temp'][i,:].squeeze()
    
    fid2 = nc.Dataset(nc2)
    sst2 = fid2.variables['temp'][i,:].squeeze()
    sst = sst2-sst1 
    print np.min(sst),np.max(sst)
    return sst,mon

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# CCS grid shape ONLY
GRD = pyroms.grid.get_ROMS_grid('CCS')
mask = GRD.hgrid.mask_rho[:]
glat = GRD.hgrid.lat_rho[:]
glon = GRD.hgrid.lon_rho[:]
dep  = GRD.vgrid.h[:]
#glat = GRD.hgrid.lat_rho[1:-1,1:-1]
#glon = GRD.hgrid.lon_rho[1:-1,1:-1]
#angs = GRD.hgrid.angle_rho[1:-1,1:-1]
#mask = GRD.hgrid.mask_rho[1:-1,1:-1]
#glat = GRD.hgrid.lat_rho[1:-1,1:-1]
#glon = GRD.hgrid.lon_rho[1:-1,1:-1]
#angs = GRD.hgrid.angle_rho[1:-1,1:-1]


### OFFSETS
joffset = 0
ioffset = 0

m_offset = 0.05
mask_val = 0
map_order = 30
vip_eta = [0,870]
vip_xi  = [0,376]

# INITIAL FIGURE
#fig, ax = plt.subplots(figsize=(4,6))
#fig, ax = plt.subplots(figsize=(4,8))
#fig, ax = plt.subplots(figsize=(4,4))
fig, ax = plt.subplots(figsize=(4,6))
fig.tight_layout()
m = Basemap(llcrnrlat=np.min(glat)-m_offset,urcrnrlat = np.max(glat)+m_offset,llcrnrlon=np.min(glon)-m_offset,urcrnrlon=np.max(glon)+m_offset, resolution='i', ax=ax)

P = m.pcolormesh(glon,glat,mask,vmin=.5,vmax=.75,edgecolors='face',cmap='Blues',zorder=map_order)
P.cmap.set_under('white')
P.cmap.set_over([.9,.97,1])

outline_mask(m,mask,mask_val,glon[0,0],glat[0,0],glon[-1,-1],glat[-1,-1])

#DOMAIN OUTLINE
for j in range(glat.shape[0]-2):
    m.plot((glon[j,0],glon[j+1,0]),(glat[j,0],glat[j+1,0]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((glon[j,-1],glon[j+1,-1]),(glat[j,-1],glat[j+1,-1]),linewidth=2,color='k',zorder=map_order+1)
for ii in range(glat.shape[1]-2):
    m.plot((glon[0,ii],glon[0,ii+1]),(glat[0,ii],glat[0,ii+1]),linewidth=2,color='k',zorder=map_order+1)
    m.plot((glon[-1,ii],glon[-1,ii+1]),(glat[-1,ii],glat[-1,ii+1]),linewidth=2,color='k',zorder=map_order+1)

polygon_patch(m,ax)

tx =plt.text(-111.5,31,'', fontsize=20,zorder=map_order+5)

# TOP PLOT
#m.drawmeridians([-128,-124], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
#m.drawparallels([42,46,50], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)

# TOP MIDDLE PLOT
# m.drawmeridians([-128,-124], labels=[0,0,0,0], linewidth=2,fmt='%d', fontsize=18,zorder=map_order+5)
# m.drawparallels([36,40], labels=[0,0,0,0], linewidth=2,fmt='%d', fontsize=18,zorder=map_order+5)

# BOTTOM MIDDLE PLOT
# m.drawmeridians([-120,-116], labels=[0,0,0,0], linewidth=2.5,fmt='%d', fontsize=18,zorder=map_order+5)
# m.drawparallels([32,36], labels=[0,0,0,0], linewidth=2.5,fmt='%d', fontsize=18,zorder=map_order+5)

# BOTTOM PLOT
m.drawmeridians([-116,-112], labels=[0,0,0,0], linewidth=1,fmt='%d', fontsize=18,zorder=map_order+5)
m.drawparallels([24,28,32], labels=[0,0,0,0], linewidth=1,fmt='%d', fontsize=18,zorder=map_order+5)

sst,mon = get_sst(0)
im1 = m.pcolor(glon,glat,sst[:],vmin=0,vmax=4,cmap='nipy_spectral',zorder=map_order)
im2 = m.contour(glon,glat,dep,[100,200],linewidths=1,linestyles='dotted',colors='white',zorder=map_order+1)
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
    im1   = m.pcolormesh(glon,glat,sst,vmin=0,vmax=4,cmap='nipy_spectral',zorder=map_order)
    im2   = m.contour(glon,glat,dep,[100,200],linewidths=.5,linestyles='dotted',colors='white',zorder=map_order+1)
    polygon_patch(m,ax)
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    tx_str = mon 
    tx.set_text(tx_str)
    # ADD Colorbar on first iteration 
    #if i == 0:
    #   cbar = m.colorbar(im1, location='bottom',size="5%", pad="3%",ticks=[])

# TOP PLOT
# ax.set_ylim(40,np.max(glat)+m_offset-.4)
# ax.set_xlim(-128.8,-123.8)

# TOP MIDDLE PLOT
#ax.set_ylim(35.5,41)
#ax.set_xlim(-125,-121.5)

# BOTTOM MIDDLE PLOT
# ax.set_ylim(31.5,35)
# ax.set_xlim(-121,-116.5)

# BOTTOM PLOT
ax.set_ylim(22.7,32)
ax.set_xlim(-117,np.max(glon)+m_offset)

ani = animation.FuncAnimation(fig, updatefig,frames=12, blit=False)
# ani.save('CCS_ROMS_SST_DIFS.gif', writer = 'ImageMagickWriter',fps=1)
# ani.save('CCS_ROMS_SST_DIFS.gif', writer = 'imagemagick',fps=1)
ani.save('CCS_ROMS_BOT_DIFS.gif', writer = 'imagemagick',fps=1)
#ani.save('CCS_ROMS_VEL_CLIM.gif', writer = 'imagemagick',fps=1)
plt.show()
