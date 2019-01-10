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
rlat = fid.variables['lat_rho'][:]
rlon = fid.variables['lon_rho'][:]
plat = fid.variables['lat_psi'][:]
plon = fid.variables['lon_psi'][:]

map_order = 30
m_offset = 0.01
mask_val = 0

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
    
    m = Basemap(llcrnrlat=np.min(plat)-m_offset,urcrnrlat = np.max(plat)+m_offset,llcrnrlon=np.min(plon)-m_offset,urcrnrlon=np.max(plon)+m_offset, resolution='i', ax=ax)
    # NOTE: pcolormesh is INCORRECT - using coordinates for value rather than corners
    P = m.pcolormesh(lon_0,lat_0,out,edgecolors='face',vmin=0,vmax=10,zorder=map_order)

    outline_mask(m,mask_rho[1:-1,1:-1],mask_val,plon[0,0],plat[0,0],plon[-1,-1],plat[-1,-1])

    plt.colorbar(P)
    plt.show()


  




