import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
 
fid1 = nc.Dataset('/Users/elizabethdrenkard/ANALYSES/CCS/CCS_scripts/ERAinterim/80x80_regrid/regridded_drowned_ERAi_msl_1981-2010_monthly_clim.nc')
fid2 = nc.Dataset('/Users/elizabethdrenkard/ANALYSES/CCS/CCS_scripts/ERAinterim/70x80_regrid/regridded_drowned_ERAi_msl_1981-2010_monthly_clim.nc') 

pair1 = fid1.variables['Pair'][:].squeeze()
pair2 = fid2.variables['Pair'][:].squeeze() 

pdif = pair2-pair1

print set(pdif.flatten())

plt.pcolor(pdif)
plt.colorbar()
plt.show()


