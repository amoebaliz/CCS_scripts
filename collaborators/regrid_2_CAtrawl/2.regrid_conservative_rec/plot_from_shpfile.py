import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import shapefile as shp 
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


def plot_to_caltrawl(mon):
    im_all=[]
    # LOOP OVER EACH SHAPE IN SHAPEFILE
    for ns in range(nshapes):
        s = sf.shape(ns)
        bbox = [coord for coord in s.bbox]
       
        # EXTRACT CORNER COORDINATES
        shp_lons = np.array((bbox[0], bbox[2]))
        shp_lats = np.array((bbox[3], bbox[1]))

        # IF INSIDE MAP SUBDOMAIN
        if ((lon_bnds[plot_num,0]<shp_lons[0]<lon_bnds[plot_num,1])  or  \
            (lon_bnds[plot_num,0]<shp_lons[1]<lon_bnds[plot_num,1])) and \
           ((lat_bnds[plot_num,0]<shp_lats[0]<lat_bnds[plot_num,1])  or  \
            (lat_bnds[plot_num,0]<shp_lats[1]<lat_bnds[plot_num,1])):

           # GET THE CLIMATOLOGICAL VALUE
           c_val = sf.records()[ns][mon+1] 

           # plot each pcolor grid cell
           im_p =  m.pcolormesh(shp_lons,shp_lats,c_val*np.ones((2,2)), \
                                vmin=1.5,vmax=4.5,cmap='nipy_spectral',     \
                                edgecolors='k',zorder=map_order) 
           # add to list of images
           im_all.append(im_p)

    return im_all
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CALIFORNIA TRAWL INFORMATION
shp_fil = 'caltrawl_GCS_Climatological_CESM_dSST'
sf = shp.Reader(shp_fil)
nshapes = len(sf.shapes())

# ROMS Grid information
grdfile = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid = nc.Dataset(grdfile)
mask_rho = fid.variables['mask_rho'][:]
rlat = fid.variables['lat_rho'][:]
rlon = fid.variables['lon_rho'][:]
plat = fid.variables['lat_psi'][:]
plon = fid.variables['lon_psi'][:]

### OFFSETS

m_offset = 0.05
mask_val = 0
map_order = 30

# CHOOSE GRID ZOOM (SOUTH to NORTH)
plot_num=2

lat_bnds = np.array(((32.35,34.49),  \
                     (34.51,37.99),\
                     (38.01,42)))

lon_bnds = np.array(((-122.05,-117.065),  \
                     (-124.4,-120.455),\
                     (-125.7,-122.8)))

fig, ax = plt.subplots(figsize=(2*abs(np.diff(lon_bnds[plot_num,:])),\
                                2*abs(np.diff(lat_bnds[plot_num,:]))))

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

mon_list =['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'] 
tx =plt.text(np.max(lon_bnds[plot_num,:])-.5,\
             np.max(lat_bnds[plot_num,:])-.25,\
             '', fontsize=20,zorder=map_order+5)

ax.set_ylim(lat_bnds[plot_num,0]-m_offset,lat_bnds[plot_num,1]+m_offset)
ax.set_xlim(lon_bnds[plot_num,0],lon_bnds[plot_num,1])

m.drawmeridians(np.arange(-125,-119+1,2), labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
m.drawparallels(np.arange(34,42+1,2), labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)

ax.set_ylim(lat_bnds[plot_num,0]-m_offset,lat_bnds[plot_num,1]+m_offset)
ax.set_xlim(lon_bnds[plot_num,0]-m_offset,lon_bnds[plot_num,1]+m_offset)

# ANIMATION
def updatefig(i):
    global im_all,tx
    print i
    # REMOVE images after first step
    if i > 0:
       for im in im_all:
           im.remove()

    im_all=plot_to_caltrawl(i)
    polygon_patch(m,ax)
 
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    tx_str = mon_list[i] 
    tx.set_text(tx_str)

ani = animation.FuncAnimation(fig, updatefig,frames=12, blit=False)
ani.save('CESM_ROMS_SST_dif_CAtrawl_regrid.gif', writer = 'imagemagick',fps=1)
