import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# Grid File
grd_fil = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid = nc.Dataset(grd_fil)

# Get the variables needed for terrain-following algorithms
h = fid.variables['h'][:] 
hc = fid.variables['hc'][:]
N = 50
s_rho = fid.variables['s_rho'][:]
Cs_r = fid.variables['Cs_r'][:]
Vtrans = 4

# LOAD zeta if available otherwise assume zeta = 0
zeta = np.zeros(h.shape)

# Allocate 3D depth array 
z_r = np.zeros((N,h.shape[0],h.shape[1]))

#iterate over each depth 'k'
for k in range(N):
    z0 = (hc * s_rho[k] + h*Cs_r[k]) / \
                          (hc + h)
    # take zeta into account
    z_r[k,:] = zeta[:] + (zeta[:] + h) * z0

# Checking layer depths: 
#	surface: layer '-1' or '49' (these are the same)
#	bottom:  layer '0' 

# depths are negative meter values and report the depth
# at the center of the layer

lat = fid.variables['lat_rho'][:]
lon = fid.variables['lon_rho'][:]
plt.pcolormesh(lon, lat, -1*z_r[-1,:].squeeze())
plt.colorbar()
plt.show() 
