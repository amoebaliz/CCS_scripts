import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import ESMF
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

def get_var(i):
    MONS = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    mon = MONS[i]

    # CONSTRUCT FILE NAME STRINGS 
    nc1 = his_dir + 'CCS-LD.HCo02Y_clim_' + str(i+1).zfill(2) + '.nc'
    nc2 = fut_dir + 'CCS-LD.FCo017_clim_' + str(i+1).zfill(2) + '.nc' 

    # IF MONTHLY FILES, NO NEED TO SPECIFY TIME INDEX
    fid1 = nc.Dataset(nc1)
    var1 = fid1.variables[ocean_var][:,dep_lev,:].squeeze()
    
    fid2 = nc.Dataset(nc2)
    var2 = fid2.variables[ocean_var][:,dep_lev,:].squeeze()
    var_dif = var2-var1 
    print np.min(var_dif),np.max(var_dif)
    return var_dif, mon

def interp_to_caltrawl(var):
    im_all=[]
    sourcefield.data[...] = var.squeeze()
    # LOOP OVER EACH SHAPE IN SHAPEFILE
    for ns in range(nshapes):
    #for ns in range(20):
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

           # DEFINE CENTRAL POINT
           c_lon = np.mean(shp_lons)
           c_lat = np.mean(shp_lats)
           # DEFINE SURROUNDING (FAKE) GRID POINTS
           x = [2*bbox[0]-c_lon, c_lon, 2*bbox[2]-c_lon]
           y = [2*bbox[1]-c_lat, c_lat, 2*bbox[3]-c_lat]

           Xn, Yn = np.meshgrid(x, y)
           dest_lon[...] = Xn
           dest_lat[...] = Yn
        
           # DEFINE INTERPOLATION FUNCTION
           destfield = ESMF.Field(destgrid, name = 'CAtrawl_delta')

           regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.BILINEAR,
                    unmapped_action = ESMF.UnmappedAction.IGNORE)
 
           # INTERPOLATE TO SINGLE BOX 
           destfield = regrid(sourcefield, destfield)

           # EXTRACT CENTER VALUE
           c_val = destfield.data[1,1] 
           #print ns, c_val
           if c_val>100:
              c_val=np.nan
           # plot each pcolor grid cell
           im_p =  m.pcolormesh(shp_lons,shp_lats,c_val*np.ones((2,2)),vmin=1.5,vmax=4.5,cmap='nipy_spectral',edgecolors='k',zorder=map_order) 
           # add to list of images
           im_all.append(im_p)

    return im_all
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CALIFORNIA TRAWL INFORMATION
shp_fil = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/collaborators/Farrah_shp/caltrawl_GCS'
sf = shp.Reader(shp_fil)
nshapes = len(sf.shapes())

# CESM Grid information
grdfile = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid = nc.Dataset(grdfile)
mask_rho = fid.variables['mask_rho'][:]
rlat = fid.variables['lat_rho'][:]
rlon = fid.variables['lon_rho'][:]
plat = fid.variables['lat_psi'][:]
plon = fid.variables['lon_psi'][:]

# PREPARE SOURCE GRID OBJECT
sourcegrid = ESMF.Grid(np.array(mask_rho.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)
source_lon = sourcegrid.get_coords(0)
source_lat = sourcegrid.get_coords(1)
source_lon[...] = rlon
source_lat[...] = rlat
sourcefield = ESMF.Field(sourcegrid, name = 'ROMS_Delta')

# PREPARE DESTINATION GRID OBJECT
destgrid = ESMF.Grid(np.array((3,3)), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)
dest_lon = destgrid.get_coords(0)
dest_lat = destgrid.get_coords(1)
#destfield = ESMF.Field(destgrid, name = 'CAtrawl_delta')

# DEFINE ROMS OUTPUT INFORMATION
his_dir='/Volumes/Abalone/CCS/his2/climatologies/'
fut_dir='/Volumes/Abalone/CCS/fut_017/climatologies/'

ocean_var='temp'
dep_lev=49 # Surface=49, Bottom=0

### OFFSETS

m_offset = 0.05
mask_val = 0
map_order = 30

# CHOOSE GRID ZOOM (SOUTH to NORTH)
plot_num=1

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

    var, mon = get_var(i)
    im_all=interp_to_caltrawl(var)
    polygon_patch(m,ax)
 
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    tx_str = mon 
    tx.set_text(tx_str)
    # ADD Colorbar on first iteration 
    #if i == 0:
    #   cbar = m.colorbar(im1, location='bottom',size="5%", pad="3%",ticks=[])

ani = animation.FuncAnimation(fig, updatefig,frames=12, blit=False)
ani.save('CCS_ROMS_SST_dif_CAtrawl_regrid.gif', writer = 'imagemagick',fps=1)
#plt.show()
