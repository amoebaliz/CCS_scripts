import ESMF
import numpy as np
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection

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

    sourcefield.data[...] = fid.variables['mask_rho'][:].squeeze()
    regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.BILINEAR,
             unmapped_action = ESMF.UnmappedAction.IGNORE)
    destfield = regrid(sourcefield, destfield)
    
    return_val = destfield.data

    # CONVERTING TO 1/0 boolean for mask_where: TRUE (1) where you want to mask
    return_val[return_val<.5]=0    # FALSE (initially land and ocean outside ROMS domain) 
    return_val[return_val>.5]=1    # TRUE  (initially water inside ROMS domain)

    # CHECK NUMBER OF POINTS USED vs. TOTAL CESM SUB-DOMAIN 
    # print 'ROMS OCEAN:', np.sum(return_val)
    # print 'CESM SIZE:', return_val.shape[0]*return_val.shape[1]

    return_val = -1*return_val + 1 # REVERSE MASK VALUES FOR where_mask 
    return_val2 = return_val.astype(np.bool) # MAKE BOOLEAN

    return return_val2


def plot_annual_maps():
    str_val = ['a','b','c']

    # Rendering
    nplot=3
    fig = plt.figure(figsize=(11,5))
    gs = gridspec.GridSpec(1,nplot)
    gs.update(wspace=0.025)

    for n in range(nplot):
        ax = plt.subplot(gs[n])
        ax.text(m_lon[0]-2, m_lat[-1]+.8, str_val[n], fontsize=14)

    plt_str ='CEMS_LE_CCS_traj'
    plt.savefig(plt_str)
       

 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# NCAR directory for Large Ensemble files
dir = '/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/'

# File name strings
his_fil_bas = 'b.e11.B20TRC5CNBDRD.f09_g16.'
fut_fil_bas = 'b.e11.BRCP85C5CNBDRD.f09_g16.'

fil_bas2 = '.pop.h.SST.'

# CCS DOMAIN
a = 245
b = 317
c = 225
d = 260

# Number of Years in timeseries
nyr = 2100-1950+1

# Number of NCAR LE endmembers
end_mem = 35

# Initializing matrix for multi-endmember timeseries
var_stor=np.zeros((ny*12,end_mem))

# EACH MONTH (if climatology)
# for nm in range(nmons):
#    print nm

    for ne in range(1,end_mem+1):
        #######################
        ## HISTORICAL VALUES ##
        #######################
        if ne == 1: 
           his_str_yr = 1850
        else:
           his_str_yr = 1920

        # File name
        his_ncfil  = dir + his_fil_bas + str(ne).zfill(3) + fil_bas2 + str(his_str_yr) + '01-' + str(200512)+ '.nc'
        fid_his = nc.Dataset(his_ncfil)

        # Total number of yrs in historical file
        nyr_h = 2005-his_str_yr+1

        # From beginning of 1950 to end of length of timeseries
        Ifyr = 12*(nyr_h-(2005-1950+1))

        # Fetch Variables
        if VAR == 'PH':
           VAR_ts = fid_his.variables[VAR][Ifyr:,a:b,c:d].squeeze()
        else:
           VAR_ts = fid_his.variables[VAR][Ifyr:,:,a:b,c:d].squeeze()
        
        ################### 
        ## FUTURE VALUES ##
        ###################
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
       #plt.pcolor(ROMS_MASK); plt.colorbar();plt.show()

    # OLD WAY OF Generating CESM mask (land + ROMS)
    # WHEN SST_F did not preserve mask: 
    # SST_H[ROMS_MASK] = np.ma.masked
    # LENS_MASK = SST_H.mask
 
    # CESM MASK is preserved by current method
    # Apply ROMS mask and spatially average
    sst_all_stor[:,ROMS_MASK] = np.ma.masked
    avg_fut_sst = np.ma.mean(np.ma.mean(sst_all_stor,axis=2),axis=1).squeeze()

    dif_all_stor[:,ROMS_MASK] = np.ma.masked
    avg_dif_sst[nm,:] = np.ma.mean(np.ma.mean(dif_all_stor,axis=2),axis=1).squeeze()
    #print avg_dif_sst

    # Generating and masking 2D fields for plotting
    sst_all_ind = np.ma.array(np.argsort(sst_all_stor,axis=0)[-1]+1)
    sst_all_ind[SST_H.mask + ROMS_MASK] = np.ma.masked

    dif_all_ind = np.ma.array(np.argsort(dif_all_stor,axis=0)+1)
    # Note: The land mask needs to be reapplied bc argsort does not retain ???_all_stor mask 
    dif_all_ind[:,(SST_H.mask + ROMS_MASK)] = np.ma.masked

    # Figure rendering
    plot_time_series()
 
    np.save('LENS_clim_difs.npy', avg_dif_sst) 

plt.show()
