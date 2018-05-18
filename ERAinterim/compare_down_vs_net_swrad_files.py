import numpy as np
import pyroms
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

map_order = 30
m_offset = 0.01
mask_val = 0

vip_eta = [0,870]
vip_xi  = [0,376]

old_fils = ['drowned_ERAi_radsw_1981-2010_monthly_clim.nc.old']
new_fils = ['drowned_ERAi_radsw_1981-2010_monthly_clim.nc']

fvars = ['swrad'] 

fid_0 = nc.Dataset(old_fils[0])
fid_1 = nc.Dataset(new_fils[0])

# ORIGINAL FILES: on the sub CCMP domain
var_0 = fid_0.variables[fvars[0]][:]
lat_0 = fid_0.variables['lat'][:]
lon_0 = fid_0.variables['lon'][:]

# NEW FILES: on the ERAi global domain
var_1 = fid_1.variables[fvars[0]][:]

for nt in range(var_0.shape[0]):

    out = var_1[nt,:].squeeze() - var_0[nt,:].squeeze()

    fig = plt.figure(figsize=(8,8))
    fig.subplots_adjust(left=.1, right=.9, bottom=0, top=1)
    ax = fig.add_subplot(111, aspect='equal', autoscale_on=False)

    m = Basemap(llcrnrlat=np.min(lat)-m_offset,urcrnrlat = np.max(lat)+m_offset,llcrnrlon=np.min(lon)-m_offset,urcrnrlon=np.max(lon)+m_offset, resolution='f', ax=ax)

    P = m.pcolormesh(lon_0,lat_0,out,edgecolors='face',vmin=0,vmax=10,zorder=map_order)

    outline_mask(m,mask,mask_val,lon[0,0],lat[0,0],lon[-1,-1],lat[-1,-1])

    plt.colorbar(P)
    plt.show()


  




