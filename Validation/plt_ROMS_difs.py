import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

fid1 = nc.Dataset('CCS_HCo02Y_SST_clim.nc')
fid2 = nc.Dataset('CCS_HCo03Y_SST_clim.nc')

sst1 = fid1.variables['temp'][:].squeeze()
sst2 = fid2.variables['temp'][:].squeeze()

print sst1.shape, sst2.shape

plt.pcolor(np.mean(sst2-sst1,axis=0),vmin=-1,vmax=1,cmap='jet')
plt.colorbar()

plt.show()

