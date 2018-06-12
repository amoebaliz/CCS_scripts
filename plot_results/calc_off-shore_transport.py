import pyroms
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# CCS GRID
GRD = pyroms.grid.get_ROMS_grid('CCS')
mask = GRD.hgrid.mask_rho[:]
lat = GRD.hgrid.lat_rho[:]
lon = GRD.hgrid.lon_rho[:]

# BOTTOM LINE TO TOP
xi_rhos  = np.array((260,235,186,282))
eta_rhos = np.array(((20,190), (200,370), (410,580), (660,830)))

pos_trans_dif = np.zeros((4,12))
neg_trans_dif = np.zeros((4,12))

for mon in range(12):
    ncfil_p = '/Users/elizabethdrenkard/Desktop/ECCWO_FILES/CCS-LD.HCo02Y_5yr_his_clim_uv.nc'
    ncfil_f = '/Users/elizabethdrenkard/Desktop/ECCWO_FILES/CCS-LD.FCo017_5yr_fut_clim_uv.nc'

    #ncfil_p = '/Volumes/Abalone/CCS/his/clim/CCS_' + str(mon+1).zfill(2) + '-clim.nc'
    #ncfil_f = '/Volumes/Abalone/CCS/016/clim/CCS_' + str(mon+1).zfill(2) + '-clim.nc'
    fid_p = nc.Dataset(ncfil_p)
    fid_f = nc.Dataset(ncfil_f)

    u_p = fid_p.variables['u'][mon,:].squeeze()
    u_f = fid_f.variables['u'][mon,:].squeeze()
    v_p = fid_p.variables['v'][mon,:].squeeze()
    v_f = fid_f.variables['v'][mon,:].squeeze()

    pos_his_u = u_p.copy()
    pos_his_u[u_p<0] = 0
    neg_his_u = u_p.copy()
    neg_his_u[u_p>0] = 0        


    pos_fut_u = u_f.copy()
    pos_fut_u[u_f<0] = 0
    neg_fut_u = u_f.copy()
    neg_fut_u[u_f>0] = 0
    for tsect in range(4):        
        print mon, tsect
        his_pos_trans, tv= pyroms.tools.section_transport_z(pos_his_u, v_p , GRD, xi_rhos[tsect], xi_rhos[tsect],eta_rhos[tsect,0], eta_rhos[tsect,1], h1=0, h2=200)    
        print 'PAST', his_pos_trans

        his_neg_trans, tv = pyroms.tools.section_transport_z(neg_his_u, v_p , GRD, xi_rhos[tsect], xi_rhos[tsect], eta_rhos[tsect,0], eta_rhos[tsect,1], h1=0, h2=200) 
        print 'PAST', his_neg_trans

        fut_pos_trans, tv = pyroms.tools.section_transport_z(pos_fut_u, v_f , GRD, xi_rhos[tsect], xi_rhos[tsect], eta_rhos[tsect,0], eta_rhos[tsect,1], h1=0, h2=200)
        print 'FUTURE', fut_pos_trans

        fut_neg_trans,tv = pyroms.tools.section_transport_z(neg_fut_u, v_f , GRD, xi_rhos[tsect], xi_rhos[tsect], eta_rhos[tsect,0], eta_rhos[tsect,1], h1=0, h2=200)
        print 'FUTURE', fut_neg_trans,tv
        print xi_rhos[tsect], xi_rhos[tsect], eta_rhos    [tsect,0], eta_rhos[tsect,1]

        pos_trans_dif[tsect,mon]=(fut_pos_trans-his_pos_trans)/(1000000)
        neg_trans_dif[tsect,mon]=-1*(fut_neg_trans-his_neg_trans)/(1000000)

fig, ax = plt.subplots(4, sharex=True, figsize=(4,6))
fig.subplots_adjust(hspace=.2)
for nt in range(4):
        ax[nt].plot([-1,13],[0,0],'-.k')
        #if nt == 0:
        #   ax[3-nt].set_ylim(-2,7)
        #else:
        #   ax[3-nt].set_ylim(-2,3)

        ax[3-nt].set_xlim(-.2,11.2) 
        ax[3-nt].plot(range(12),  pos_trans_dif[nt,:],'-o', mec='k',linewidth=3,markersize=5)
        ax[3-nt].plot(range(12),  neg_trans_dif[nt,:],'-o', mec='k',linewidth=3,markersize=5)
        ax[3-nt].set_xticks([1,4,7,10]) 
        ax[3-nt].set_xticklabels([])

ax[0].set_ylim(-1,1.1)
ax[1].set_ylim(-1.2,2.5)
ax[2].set_ylim(-1.9,2.2)
ax[3].set_ylim(-2.5,4.5)

plt.savefig('TRANSPORTS_OFF')
plt.show()
