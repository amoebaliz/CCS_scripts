import numpy as np
import netCDF4 as nc
import pyroms
import ESMF
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection


def polygon_patch(mapid,axs):
    mapid.drawcoastlines(linewidth=0,zorder=36)
    mapid.drawmapboundary(fill_color=[.9,.97,1])
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor=(1,1,1), closed=False)
    axs.add_collection(lc)


def define_regrid():
    # AVHRR GRID
    lon[lon>180]=lon[lon>180]-360
    Xn, Yn = np.meshgrid(lon,lat) # AVHRR = 1D lat/lon

    destgrid = ESMF.Grid(np.array(Xn.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

    dest_lon = destgrid.get_coords(0)
    dest_lat = destgrid.get_coords(1)

    dest_lon[...] = Xn
    dest_lat[...] = Yn

    destfield = ESMF.Field(destgrid, name = 'AVHRR Sub')

    regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.BILINEAR,
             unmapped_action = ESMF.UnmappedAction.IGNORE)

    return regrid, destfield

def make_mask(regrid,destfield):
    # MASK ROMS REGION
    sourcefield.data[...] = r_mask

    avh_mask = regrid(sourcefield, destfield).data
    # convert to 1s and 0s
    avh_mask[avh_mask<.5]=0
    avh_mask[avh_mask>.5]=1
    # reverse masking convention
    avh_mask = -1*avh_mask + 1
    # make boolean
    avh_mask = avh_mask.astype(np.bool)

    return avh_mask

def plot_annual_maps():
    # Map specs
    m_lat = [20,50]; m_lon = [-140,-110]; m_off=3
    str_val = ['a','b']
    # Customize colormap
    #disc_cmap = cmap_discretize(matplotlib.cm.nipy_spectral, end_mem)
    map_order=30

    # Rendering
    nplot=2
    fig = plt.figure(figsize=(11,5))
    gs = gridspec.GridSpec(1,nplot)
    gs.update(wspace=0.025)

    for n in range(nplot):
        ax = plt.subplot(gs[n])
        m = Basemap(llcrnrlat=m_lat[0]-m_off,urcrnrlat=m_lat[1]+m_off,\
            llcrnrlon=m_lon[0]-m_off,urcrnrlon=m_lon[1]+m_off, resolution='i', ax=ax)
        P = m.pcolormesh(lon,lat,stor_vars[n,:].squeeze(),vmin=-3,vmax=3,cmap='jet')
        polygon_patch(m,ax)
        m.drawmeridians(m_lon, labels=[0,0,1,0], fmt='%d', fontsize=14)
        m.drawparallels(m_lat, labels=[0,0,0,0], fmt='%d')
        ax.text(m_lon[0]-2, m_lat[-1]+.8, str_val[n], fontsize=14)
        if n:
           # plot colorbar axis
           bbox = ax.get_position()
           cax = fig.add_axes([bbox.xmax*1.0+.02, bbox.ymin, bbox.width*0.08, bbox.height-.085])
           cb = plt.colorbar(P,cax=cax,ticks=np.linspace(-3,3,3))
           cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=12)
           cb.ax.set_title('$^\circ$C',horizontalalignment='center',\
                           multialignment='center')
           #plt_str = 'SODA_vs_AVHRR:avg_diff_RMS'
           plt_str = 'CCS_ROMS_vs_AVHRR:avg_diff_RMS'
           plt.savefig(plt_str)
           plt.show()
        else:
           m.drawparallels(m_lat, labels=[1,0,0,0], fmt='%d',fontsize=14)


############
# ROMS GRID
GRD = pyroms.grid.get_ROMS_grid('CCS')
r_mask = GRD.hgrid.mask_rho
r_lat = GRD.hgrid.lat_rho
r_lon = GRD.hgrid.lon_rho

sourcegrid = ESMF.Grid(np.array(r_lon.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

source_lon = sourcegrid.get_coords(0)
source_lat = sourcegrid.get_coords(1)

source_lon[...] = r_lon
source_lat[...] = r_lat

sourcefield = ESMF.Field(sourcegrid, name = 'CCS ROMS')

###########
# SODA BIAS
#soda_bias = -2.75704143254

# ACCESS ROMS SST CLIM
#SST_roms = nc.Dataset('/Users/elizabethdrenkard/TOOLS/CCS_scripts/Validation/CCS_HCo02Y_SST_clim.nc').variables['temp'][:].squeeze()
#SST_roms = SST_roms - soda_bias
 
# SODA SST CLIM
SST_roms = nc.Dataset('soda_clim_sst.nc').variables['temp'][:].squeeze()
avh_dir = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/Validation/'

#CCMP YEAR
#avh_dir = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/Validation/CCMP_check/'

for nmon in range(12):
    avh_fil = avh_dir + 'AVHRR_SST_clim_' + str(nmon+1).zfill(2) + '_clim.nc' 
    #avh_fil = avh_dir + 'SST_all_' + str(nmon+1).zfill(2) + '.nc'
    fid_avh = nc.Dataset(avh_fil)
    avh_sst = fid_avh.variables['BSST'][:].squeeze()
    
    if nmon==0:
       lat = fid_avh.variables['lat'][:]
       lon = fid_avh.variables['lon'][:]
       regrid, destfield = define_regrid()
       avh_mask = make_mask(regrid,destfield)
       squar_err = np.ma.zeros((12,avh_mask.shape[0],avh_mask.shape[1])) 
       diff_vals = np.ma.zeros((12,avh_mask.shape[0],avh_mask.shape[1]))
       stor_vars = np.ma.zeros((2,avh_mask.shape[0],avh_mask.shape[1])) 
    # calculate 'squared error'
    sourcefield.data[...] = SST_roms[nmon,:].squeeze()
    diff_vals[nmon,:] = (regrid(sourcefield, destfield).data - avh_sst)
    squar_err[nmon,:] = (regrid(sourcefield, destfield).data - avh_sst)**2

    #diff_vals[nmon,avh_mask]=np.ma.masked
    #plt.pcolor(lon,lat,diff_vals[nmon,:].squeeze(),vmax=1,cmap='jet'); plt.colorbar(); plt.show()

rms_err = np.sqrt(np.mean(squar_err,axis=0))
avg_dif = np.ma.mean(diff_vals,axis=0)

stor_vars[0,:]=avg_dif
stor_vars[1,:]=rms_err
stor_vars[:,avh_mask] = np.ma.masked

# DETERMINING AVERAGE SODA BIAS IN CCS DOMAIN
#avg_dif[avh_mask] = np.ma.masked
#avg_dif[np.abs(avg_dif)>100] = np.ma.masked
#print np.ma.mean(avg_dif)
plot_annual_maps()


