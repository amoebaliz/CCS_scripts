import subprocess
import os
import sys 

import pyroms
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
        if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
           v.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))
        else :
           l.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))

        if p[1] == mask_img.shape[1] - 1 :
           if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
               v.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))
           else:
               l.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))
        else :
           if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
              v.append((lon[p[0]+1,p[1]+1],lat[p[0]+1,p[1]+1]))
           else:
              l.append((lon[p[0]+1,p[1]+1],lat[p[0]+1,p[1]+1]))

        l.append((np.nan,np.nan))
        v.append((np.nan,np.nan))
    #vertical segments
    for p in zip(*ver_seg):
        if p[1] == mask_img.shape[1]-1:
           if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
              v.append((lon[p[0],p[1]],lat[p[0],p[1]]))
              v.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))
           else:
              l.append((lon[p[0],p[1]],lat[p[0],p[1]]))
              l.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))
        elif p[0] == mask_img.shape[0]-1:
             if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
              v.append((lon[p[0],p[1]],lat[p[0],p[1]]))
              v.append((lon[p[0]+1,p[1]],lat[p[0]+1,p[1]]))
             else:
              l.append((lon[p[0],p[1]],lat[p[0],p[1]]))
              l.append((lon[p[0],p[1]+1],lat[p[0],p[1]+1]))
        else:
           if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
              v.append((lon[p[0],p[1]+1],lat[p[0],p[1]+1]))
              v.append((lon[p[0]+1,p[1]+1],lat[p[0]+1,p[1]+1]))
           else:
              l.append((lon[p[0],p[1]+1],lat[p[0],p[1]+1]))
              l.append((lon[p[0]+1,p[1]+1],lat[p[0]+1,p[1]+1]))

        l.append((np.nan, np.nan))
        v.append((np.nan, np.nan))
    segments = np.array(l)
    vip_segments = np.array(v)
    mapid.plot(segments[:,0], segments[:,1], latlon=True, color=(0,0,0), linewidth=.75,zorder=map_order+2)
    mapid.plot(vip_segments[:,0], vip_segments[:,1], latlon=True, color=(0,0,0), linewidth=.75,zorder=map_order+3)


# CCS DETAILS

GRD = pyroms.grid.get_ROMS_grid('CCS')
mask = GRD.hgrid.mask_rho
lat = GRD.hgrid.lat_rho
lon = GRD.hgrid.lon_rho

era_dir1 = '/Users/liz.drenkard/external_data/ERAinterim/drowned/original_drowned_ERAi/'
era_dir2 = '/Users/liz.drenkard/external_data/ERAinterim/drowned/'

filelst = subprocess.check_output(['ls', era_dir2]).replace('/n',' ').split()

varlst = ['Pair', 'rain', 'Qair','lwrad_down','swrad','Tair']

map_order = 30
m_offset = 0.01
mask_val = 0

vip_eta = [0,870]
vip_xi  = [0,376]

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

        m = Basemap(llcrnrlat=np.min(lat)-m_offset,urcrnrlat = np.max(lat)+m_offset,llcrnrlon=np.min(lon)-m_offset,urcrnrlon=np.max(lon)+m_offset, resolution='f', ax=ax)
        
        P = m.pcolormesh(lon_e,lat_e,dif_val,edgecolors='face',vmin=-10, vmax=10,zorder=map_order)

        outline_mask(m,mask,mask_val,lon[0,0],lat[0,0],lon[-1,-1],lat[-1,-1])
        
        plt.colorbar(P)

    plt.show()

