import subprocess
import os
import sys 

import numpy as np
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt 

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
plat = fid.variables['lat_psi'][:]
plon = fid.variables['lon_psi'][:]

era_dir1 = '/Users/liz.drenkard/external_data/ERAinterim/drowned/original_drowned_ERAi/'
era_dir2 = '/Users/liz.drenkard/external_data/ERAinterim/drowned/'

filelst = subprocess.check_output(['ls', era_dir2]).replace('/n',' ').split()

varlst = ['Pair', 'rain', 'Qair','lwrad_down','swrad','Tair']

map_order = 30
m_offset = 0.01
mask_val = 0

for (fil,var) in zip(filelst,varlst):
    print var 
    ncfil_e1 = era_dir1 + fil
    ncfil_e2 = era_dir2 + fil

    fid_e1 = nc.Dataset(ncfil_e1)
    fid_e2 = nc.Dataset(ncfil_e2)

    var_e1 = fid_e1.variables[var][:]
    var_e2 = fid_e2.variables[var][:]
    
    lat_e = fid_e1.variables['lat'][:]
    lon_e = fid_e1.variables['lon'][:]

    for nt in range(var_e1.shape[0]):
        dif_val = var_e2[nt,:].squeeze()-var_e1[nt,:].squeeze()

        fig = plt.figure(figsize=(8,8))
        fig.subplots_adjust(left=.1, right=.9, bottom=0, top=1)
        ax = fig.add_subplot(111, aspect='equal', autoscale_on=False)

        m = Basemap(llcrnrlat=np.min(plat)-m_offset,urcrnrlat = np.max(plat)+m_offset,llcrnrlon=np.min(plon)-m_offset,urcrnrlon=np.max(plon)+m_offset, resolution='i', ax=ax)
        # NOTE: pcolormesh is INCORRECT bc using coordinates of value location rather than edges 
        P = m.pcolormesh(lon_e,lat_e,dif_val,edgecolors='face',vmin=-10, vmax=10,zorder=map_order)

        outline_mask(m,mask_rho[1:-1,1:-1],mask_val,plon[0,0],plat[0,0],plon[-1,-1],plat[-1,-1])
        
        plt.colorbar(P)

    plt.show()

