import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt


fidu = nc.Dataset('/glade/p/work/edrenkar/external_data/LENS/difs/016_U_clim_delta.nc')  
fidv = nc.Dataset('/glade/p/work/edrenkar/external_data/LENS/difs/016_V_clim_delta.nc')
fidu10 = nc.Dataset('/glade/p/work/edrenkar/external_data/LENS/scripts/difs/016_U10_clim_delta.nc')

u = fidu.variables['U'][:].squeeze()
v = fidv.variables['V'][:].squeeze()
u10 = fidu10.variables['U10'][:].squeeze()

S = np.sqrt(u**2+v**2)
for nt in range(12):
    plt.figure()
    plt.pcolor(S[nt,:,:].squeeze(),vmin=-1,vmax=1)
    plt.colorbar()

    plt.figure()
    plt.pcolor(u10[nt,:,:].squeeze(),vmin=-1,vmax=1)
    plt.colorbar()
    plt.show()

