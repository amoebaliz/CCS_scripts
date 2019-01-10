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

def get_sst(i):
    MONS = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    mon = MONS[i] 
    ncfile = '/Users/elizabethdrenkard/Desktop/soda_fils/soda3.4.1_1981-2010_clim_'+str(i+1).zfill(2) + '.nc'
    fid = nc.Dataset(ncfile)
    sst = fid.variables['temp'][0,0,:].squeeze()    
    sst = fill_CA_Gulf(sst)
    lat = fid.variables['yt_ocean'][:] 
    lon = fid.variables['xt_ocean'][:]   
    return sst,lat,lon, mon

def get_vel(i):
    MONS = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    mon = MONS[i]
    ncfile = '/Users/elizabethdrenkard/Desktop/soda_fils/soda3.4.1_1981-2010_clim_'+str(i+1).zfill(2) + '.nc'
    fid = nc.Dataset(ncfile)
    u = fid.variables['u'][0,0,:].squeeze()
    u = fill_CA_Gulf(u)
    v = fid.variables['v'][0,0,:].squeeze()
    v = fill_CA_Gulf(v)
    mag = np.sqrt(u**2+v**2)
    lat = fid.variables['yu_ocean'][:]
    lon = fid.variables['xu_ocean'][:]
    return u, v, mag, lat, lon, mon

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# DATA LOCATION 
dir = '/Users/elizabethdrenkard/Desktop/soda_fils/'
fil_pre = 'soda3.4.1_1981-2010_clim_'

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
afreq = 8
#plt.title('SST (oC)')
tx =plt.text(-118,46.5,'', fontsize=20,zorder=map_order+5)
n=0

m.drawmeridians([-142,-111], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
m.drawparallels([18,50], labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)

# NOTE: this use of pcolor is incorrect - lats/lons should be the corners, not the actual location
u,v,mag,lats,lons,mon = get_vel(0)
im1 = m.pcolor(lons,lats,mag,vmin=0,vmax=.1,cmap='Oranges',zorder=map_order)
im2 = m.quiver(lons[::afreq],lats[::afreq],u[::afreq,::afreq],v[::afreq,::afreq],zorder=map_order+2)

#plt.show()
# ANIMATION
def updatefig(i):
    global im1,im2,tx
    print i
    # REMOVE images after first step
    im1.remove()
    im2.remove()
    #if i > 0:
    #   im1.remove()
    #sst,lats,lons,mon = get_sst(i)
    u,v,mag,lats,lons,mon = get_vel(i)
    #im1   = m.pcolormesh(lons,lats,sst,vmin=8,vmax=28,cmap='nipy_spectral',zorder=map_order)
    im1 = m.pcolor(lons,lats,mag,vmin=0,vmax=.1,cmap='OrRd',zorder=map_order)
    im2 = m.quiver(lons[::afreq],lats[::afreq],u[::afreq,::afreq],v[::afreq,::afreq],scale=5,zorder=map_order+2)
    qk = m.quiverkey(im2,-118,46.5,.1,r'$.1 \frac{m}{s}$', labelpos='W',zorder=map_order+2)
    polygon_patch(m,ax)
    tx_str = mon 
    tx.set_text(tx_str)
    # ADD Colorbar on first iteration 
    if i == 0:
       cbar = m.colorbar(im1, location='bottom',size="5%", pad="3%",ticks=[])

ani = animation.FuncAnimation(fig, updatefig,frames=12, blit=False)
#ani.save('SODA_SST_CLIM.gif', writer = 'imagemagick',fps=1)
ani.save('SODA_VEL_CLIM.gif', writer = 'imagemagick',fps=1)
plt.show()
