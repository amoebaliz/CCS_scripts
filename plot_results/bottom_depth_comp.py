import pyroms
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection


# CCS grid shape ONLY
GRD = pyroms.grid.get_ROMS_grid('CCS')
roms_mask = GRD.hgrid.mask_rho[:]
glat = GRD.hgrid.lat_rho[:]
glon = GRD.hgrid.lon_rho[:]
dep  = GRD.vgrid.h[:]

MASK = np.ones(dep.shape,dtype=bool)
MASK[dep>=30]=0
MASK[roms_mask==0]=0

print MASK.shape
dep = np.ma.masked_array(dep)
dep[MASK==0]=np.ma.masked

#plt.pcolor(dep)
#plt.colorbar()
#plt.show()

# nc_files
his_fil='/glade/p/work/edrenkar/MODELS/CCS/ANALYSES/CCS-LD.HCo02Y/5yr_his_BOT_TEMP_clim.nc'
fut_fil='/glade/p/work/edrenkar/MODELS/CCS/ANALYSES/CCS-LD.FCo017/5yr_fut_BOT_TEMP_clim.nc'

h_fid = nc.Dataset(his_fil)
f_fid = nc.Dataset(fut_fil)

h_temp = h_fid.variables['temp'][:].squeeze()
f_temp = f_fid.variables['temp'][:].squeeze()

temp_dif = f_temp-h_temp

avg_dif = np.mean(temp_dif,axis=0)
avg_dif = np.ma.masked_where(MASK==0, avg_dif)

print 'MEEEEEP', np.ma.mean(avg_dif)

for nm in range(temp_dif.shape[0]):
    m_temp_dif = temp_dif[nm,:].squeeze()
    m_temp_dif = np.ma.masked_where(MASK==0, m_temp_dif)
    print np.ma.mean(m_temp_dif)
 
    m_his_temp = h_temp[nm,:].squeeze()
    m_his_temp = np.ma.masked_where(MASK==0, m_his_temp)
    print np.ma.mean(m_his_temp)

    m_fut_temp = f_temp[nm,:].squeeze()
    m_fut_temp = np.ma.masked_where(MASK==0, m_fut_temp)
    print np.ma.mean(m_fut_temp)
#    print 'STD:' ,np.ma.std(m_fut_temp) 

#temp_dif1 = np.ma.masked_where(MASK==1,temp_dif)

#print np.ma.max(np.ma.masked_where(MASK==1,temp_dif))
#print np.ma.min(np.ma.masked_where(MASK==1,temp_dif))

# print np.ma.mean(h_temp[:,MASK])#,axis=(2,1))
# print np.ma.mean(f_temp[:,MASK])#,axis=(2,1))
# print np.mean(f_temp[:,MASK]-h_temp[:,MASK])

