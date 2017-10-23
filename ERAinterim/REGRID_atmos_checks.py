import numpy as np
import netCDF4 as nc
import pyroms
import regrid_atmos
import matplotlib.pyplot as plt

grd = pyroms.grid.get_ROMS_grid('CCS')
Yout = np.asfortranarray(grd.hgrid.lat_rho.astype(float))
Xout = np.asfortranarray(grd.hgrid.lon_rho.astype(float))

#SWFSC
#fid1 = nc.Dataset('/Users/liz.drenkard/external_data/ERAinterim/drowned/80_x_80/drowned_ERAi_msl_1981-2010_monthly_clim.nc')
#fid2 = nc.Dataset('/Users/liz.drenkard/external_data/ERAinterim/drowned/70_x_70/drowned_ERAi_msl_1981-2010_monthly_clim.nc')

#MACBOOK
fid1 = nc.Dataset('/Users/elizabethdrenkard/external_data/ERAinterim/drowned/80_x_80/drowned_ERAi_msl_1981-2010_monthly_clim.nc')
fid2 = nc.Dataset('/Users/elizabethdrenkard/external_data/ERAinterim/drowned/70_x_80/drowned_ERAi_msl_1981-2010_monthly_clim.nc')

p1 = fid1.variables['Pair'][0,:].squeeze()
p2 = fid2.variables['Pair'][0,:].squeeze()

Amin1 = np.min(p1)
Amin2 = np.min(p2)

Amax1 = np.max(p1)
Amax2 = np.max(p2)

lat1 = fid1.variables['lat'][:].squeeze().astype(float)
lat2 = fid2.variables['lat'][:].squeeze().astype(float)

lon1 = fid1.variables['lon'][:].squeeze().astype(float)
lon1[lon1>180]=lon1[lon1>180]-360

lon2 = fid2.variables['lon'][:].squeeze().astype(float)
lon2[lon2>180]=lon2[lon2>180]-360

Xinp1 = np.asfortranarray(np.tile(lon1,(len(lat1),1))) 
Xinp2 = np.asfortranarray(np.tile(lon2,(len(lat2),1))) 

Yinp1 = np.asfortranarray(np.transpose(np.tile(lat1,(len(lon1),1))))
Yinp2 = np.asfortranarray(np.transpose(np.tile(lat2,(len(lon2),1))))

Jout1, Iout1, Fout1 = regrid_atmos.regrid_atmos(Xinp1, Yinp1, np.asfortranarray(p1.squeeze().astype(float)), Amin1, Amax1, Xout, Yout)
Jout2, Iout2, Fout2 = regrid_atmos.regrid_atmos(Xinp2, Yinp2, np.asfortranarray(p2.squeeze().astype(float)), Amin2, Amax2, Xout, Yout)
store_J_dif = np.zeros(Jout1.shape[1])
for nt in range(Jout1.shape[1]):
    store_J_dif[nt] = Jout1[nt,nt]-Jout2[nt,nt]-10

print np.sum(abs(store_J_dif))
print Jout1.shape
print 873*382
print np.min(Jout1-Jout2-10)
print np.max(Jout1-Jout2-10)
plt.figure()
plt.pcolor(Fout1-Fout2)
plt.colorbar()
#print np.min(Fout1-Fout2), np.max(Fout1-Fout2)

#plt.figure()
#plt.pcolor(Iout1-Iout2)
#plt.colorbar()
#print np.min(Iout1-Iout2), np.max(Iout1-Iout2)

plt.figure()
plt.pcolor(Jout1-Jout2-10)
plt.colorbar()
#print np.min(Jout1-Jout2), np.max(Jout1-Jout2)

plt.show()
