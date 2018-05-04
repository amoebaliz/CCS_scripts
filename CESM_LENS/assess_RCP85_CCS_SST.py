import ESMF
import numpy as np
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

def ccs_roms_region(lat,lon):
    # CCS ROMS GRID
    # without pyroms
    fid = nc.Dataset('/glade/p/work/edrenkar/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc')
    r_lat = fid.variables['lat_rho'][:]
    r_lon = fid.variables['lon_rho'][:]
    sourcegrid = ESMF.Grid(np.array(r_lon.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG) 

    # CESM GRID
    lon[lon>180]=lon[lon>180]-360 
    destgrid = ESMF.Grid(np.array(lon.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

    ## POINTERS
    source_lon = sourcegrid.get_coords(0)
    source_lat = sourcegrid.get_coords(1)
    dest_lon = destgrid.get_coords(0)
    dest_lat = destgrid.get_coords(1)

    ## FILLS
    source_lon[...] = r_lon
    source_lat[...] = r_lat
    dest_lon[...] = lon
    dest_lat[...] = lat

    sourcefield = ESMF.Field(sourcegrid, name = 'CCS ROMS')
    destfield = ESMF.Field(destgrid, name = 'CESM Sub')

    sourcefield.data[...] = fid.variables['h'][:].squeeze()
    regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.BILINEAR,
             unmapped_action = ESMF.UnmappedAction.IGNORE)
    destfield = regrid(sourcefield, destfield)
    
    return_val = destfield.data

    # CONVERTING TO 1/0 boolean for mask_where
    return_val[return_val==0]=1  # TRUE MASK 
    return_val[return_val>1]=0   # FALSE MASK
    return_val2 = return_val.astype(np.bool) # MAKE BOOLEAN

    return return_val2

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet. 
        N: number of colors.
    """
    if type(cmap) == str:
        cmap = get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki, key in enumerate(('red','green','blue')):
        cdict[key] = [(indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in range(N+1)]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

def plot_annual_maps():
    # Map specs
    m_lat = [20,50]; m_lon = [-140,-110]; m_off=3
    str_val = ['a','b']
    # Customize colormap
    disc_cmap = cmap_discretize(matplotlib.cm.nipy_spectral, end_mem)
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
        P = m.pcolormesh(lon,lat,dif_all_ind[-1*n,:].squeeze(),vmin=.5,vmax=end_mem+.5,cmap=disc_cmap)
        polygon_patch(m,ax)
        m.drawmeridians(m_lon, labels=[0,0,1,0], fmt='%d', fontsize=14)
        m.drawparallels(m_lat, labels=[0,0,0,0], fmt='%d')
        ax.text(m_lon[0]-2, m_lat[-1]+.8, str_val[n], fontsize=14)
        if n:
           print 'MEEP' 
           # plot colorbar axis
           bbox = ax.get_position()
           cax = fig.add_axes([bbox.xmax*1.0+.02, bbox.ymin, bbox.width*0.08, bbox.height-.085])
           cb = plt.colorbar(P,cax=cax,ticks=np.linspace(5,35,7))
           cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=12)
           cb.ax.set_title('LENS\nEndmember',horizontalalignment='center',\
                           multialignment='center')
           plt_str = str(nm).zfill(2) + '_MIN_DIF_SST_ENDMEMBER'
           plt.savefig(plt_str)
        else:
           m.drawparallels(m_lat, labels=[1,0,0,0], fmt='%d',fontsize=14)
        
def plot_clim_difs():
     # NOTE #5: Line plot
     print avg_dif_sst.shape
     # Plot each endmember against month
     #for nl in range(avg_dif_sst.shape[1]): 
     #    plt.plot(np.arange(1,12+1),avg_dif_sst[:,nl].squeeze())
     #    plt.show()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# NCAR directory for Large Ensemble files
dir = '/glade/p/cesmLE/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/SST/'

# File name strings
his_fil_bas = 'b.e11.B20TRC5CNBDRD.f09_g16.'
fut_fil_bas = 'b.e11.BRCP85C5CNBDRD.f09_g16.'

fil_bas2 = '.pop.h.SST.'

# CCS DOMAIN
a = 245
b = 317
c = 225
d = 260

# WHETHER CALCULATING CLIMATOLOGY OR ANNUAL MEAN
# 0 = annual mean; 1 = monthly climatology
CLIM=1

if CLIM:
   nmons = 12
else:
   nmons = 1
print nmons
# Number of Years in climatology
nclim_yr = 25

# Number of NCAR LE endmembers
end_mem = 35

# Customize colormap
disc_cmap = cmap_discretize(matplotlib.cm.nipy_spectral, end_mem)
map_order=30

# Initializing difference matrix for line plot
avg_dif_sst=np.zeros((nmons,end_mem))

# EACH MONTH (if climatology)
for nm in range(nmons):
    print nm

    # STORAGE ARRAYS FOR MONTHLY CLIM MEAN/STD
    sst_all_stor = np.empty((0,b-a,d-c))
    dif_all_stor = np.empty((0,b-a,d-c))

    for ne in range(1,end_mem+1):
        ############################
        ## HISTORICAL CLIMATOLOGY ##
        ############################
        if ne == 1: 
           his_str_yr = 1850
        else:
           his_str_yr = 1920

        # File name
        his_ncfil  = dir + his_fil_bas + str(ne).zfill(3) + fil_bas2 + str(his_str_yr) + '01-' + str(200512)+ '.nc'
        fid_his = nc.Dataset(his_ncfil)

        # Montly indicies for historical fields
        nyr_h = 2005-his_str_yr+1
        Ih = np.arange(nm + 12*(nyr_h-nclim_yr), 12*nyr_h, nmons) 

        # Historical climatology
        SST_H = np.mean(fid_his.variables['SST'][Ih,:,a:b,c:d].squeeze(), axis=0)
        
        ######################## 
        ## FUTURE CLIMATOLOGY ##
        ########################
        fut_str_yr = [2006,2081]
        fut_end_yr = [2080,2100]

        if ne < 34:
           # Two files: 2006-2080, 2081-2100
           tot_yr_f = [nclim_yr-np.diff(fut_end_yr),np.diff(fut_end_yr)]
        else:
           # One file: 2006-2100
           tot_yr_f = [nclim_yr]

        # Opening and appending variable number of files
        SST_F = np.empty((0,b-a,d-c))
        nfils = len(tot_yr_f)
        for nf in range(nfils):
            fut_ncfil = dir + fut_fil_bas + str(ne).zfill(3) + fil_bas2 + str(fut_str_yr[nf]) + '01-' + str(fut_end_yr[nf-nfils])+ '12.nc'
            fid_fut = nc.Dataset(fut_ncfil)
            # Number of years contained in future file 
            nyr_f = fut_end_yr[nf]-fut_str_yr[nf-nfils]+1
            # (Montly) indicies for future fields
            If = np.arange(nm + int(12*(nyr_f-tot_yr_f[nf])), 12*nyr_f, nmons)
            # Monthly climatology/Annual Mean
            SST_F = np.ma.append(SST_F, fid_fut.variables['SST'][If,:,a:b,c:d].squeeze(), axis=0)
        # Climatology
        SST_F = np.mean(SST_F,axis=0)

        # STORE VALUES FOR ANALYSIS 
        sst_all_stor = np.ma.append(sst_all_stor, np.ma.expand_dims(SST_F,axis=0), axis=0)
        dif_all_stor = np.ma.append(dif_all_stor, np.ma.expand_dims(SST_F-SST_H,axis=0), axis=0)
    ###########
    # MASKING #
    ###########
    if nm == 0:
       # Generate ROMS mask
       lat = fid_his.variables['TLAT'][a:b,c:d]
       lon = fid_his.variables['TLONG'][a:b,c:d]
       ROMS_MASK = ccs_roms_region(lat,lon)

    # Generate CESM mask (land + ROMS) 
    #SST_H[ROMS_MASK] = np.ma.masked
    #LENS_MASK = SST_H.mask

    # Apply masks and spatially average
    sst_all_stor[:,ROMS_MASK] = np.ma.masked
    avg_fut_sst = np.ma.mean(np.ma.mean(sst_all_stor,axis=2),axis=1).squeeze()

    dif_all_stor[:,ROMS_MASK] = np.ma.masked
    avg_dif_sst[nm,:] = np.ma.mean(np.ma.mean(dif_all_stor,axis=2),axis=1).squeeze()
    print avg_dif_sst
    ##########
    # OUTPUT #
    ##########
    #print "COLDEST"
    #print 'ABS', np.argsort(avg_fut_sst)[:5]+1
    #print 'DIF', np.argsort(avg_dif_sst)[:5]+1

    #print "WARMEST"
    #print 'ABS', np.argsort(avg_fut_sst)[-5:]+1
    #print 'DIF', np.argsort(avg_dif_sst)[-5:]+1

    #print "MIDDLE"
    #print np.argsort(sst_stor)[(33/2)-2:(33/2)+3]+1   
    #print np.argsort(dif_stor)[(33/2)-2:(33/2)+3]+1

    ###########
    # FIGURES #
    ###########

    # Generating and masking 2D fields for plotting
    # sst_all_ind = np.ma.array(np.argsort(sst_all_stor,axis=0)[-1]+1)
    # sst_all_ind[SST_H.mask + ROMS_MASK] = np.ma.masked

    # dif_all_ind = np.ma.array(np.argsort(dif_all_stor,axis=0)+1)
    # Note: The land mask needs to be reapplied bc argsort does not retain ???_all_stor mask 
    # dif_all_ind[:,(SST_H.mask + ROMS_MASK)] = np.ma.masked

    # Figure rendering
    # if nmons == 1:
    #  plot_annual_maps()

# LINE PLOT
if nmons>1:
   np.save('LENS_clim_difs.npy', avg_dif_sst) 
   #plot_clim_difs()


#ax.plot(np.arange(1,33+1),sst_stor,'o')
#ax.errorbar(np.arange(1,33+1),sst_stor,yerr=sst_std,fmt='o')   
#ax.plot([1,33],sst_stor[6-1]*np.ones(2),'-b')
#ax.plot([1,33],sst_stor[16-1]*np.ones(2),'-r')
#ax.plot([1,33],sst_stor[33-1]*np.ones(2),'-k')
#plt.plot(SST)
    
#ax.errorbar(np.arange(1,33+1),sst_stor,yerr=sst_std,fmt='o')
#ax.plot(np.arange(1,33+1),sst_stor,'o')
#plt.show()
