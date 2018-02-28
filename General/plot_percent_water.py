import numpy as np
import pyroms
import matplotlib.pyplot as plt

grd = pyroms.grid.get_ROMS_grid('CCS')

h = grd.vgrid.h[:]
z = grd.vgrid.z_w[:]
dz = np.diff(z,axis=0)

perc_dz = (100*dz/h).reshape((50,-1))

avg_perc = np.mean(perc_dz,axis=1)
std_perc = np.std(perc_dz,axis=1)


f, ax = plt.subplots(1,figsize=(4,8))

ax.errorbar(avg_perc, range(50),xerr=std_perc, fmt='o')
ax.plot(avg_perc,range(50),'go',zorder=3)
ax.set_xlabel('Percent')
ax.set_ylabel('Model Layer')
plt.savefig('avg_water_col_perc')
plt.show()





