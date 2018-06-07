import pyroms
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# CCS GRID
GRD = pyroms.grid.get_ROMS_grid('CCS')
mask = GRD.hgrid.mask_rho[:]
lat = GRD.hgrid.lat_rho[:]
lon = GRD.hgrid.lon_rho[:]

eta_rhos = [64,325,580,804]
xi_rhos  = np.array(((175,345),(150,320),(101,271),(197,367)))
pos_trans_dif = np.zeros((4,12))
neg_trans_dif = np.zeros((4,12))

for mon in range(12):
    #ncfil_p = '/Volumes/Abalone/CCS/his/clim/CCS_' + str(mon+1).zfill(2) + '-clim.nc'
    #ncfil_f = '/Volumes/Abalone/CCS/016/clim/CCS_' + str(mon+1).zfill(2) + '-clim.nc'

    ncfil_p = '/glade/p/work/edrenkar/MODELS/CCS/ANALYSES/CCS-LD.HCo02Y/CCS-LD.HCo02Y_5yr_his_clim_uv.nc'
    ncfil_f = '/glade/p/work/edrenkar/MODELS/CCS/ANALYSES/CCS-LD.FCo017/CCS-LD.FCo017_5yr_fut_clim_uv.nc'

    fid_p = nc.Dataset(ncfil_p)
    fid_f = nc.Dataset(ncfil_f)

    u_p = fid_p.variables['u'][mon,:].squeeze()
    u_f = fid_f.variables['u'][mon,:].squeeze()
    v_p = fid_p.variables['v'][mon,:].squeeze()
    v_f = fid_f.variables['v'][mon,:].squeeze()

    pos_his_v = v_p.copy()
    pos_his_v[v_p<0] = 0
    neg_his_v = v_p.copy()
    neg_his_v[v_p>0] = 0        


    pos_fut_v = v_f.copy()
    pos_fut_v[v_f<0] = 0
    neg_fut_v = v_f.copy()
    neg_fut_v[v_f>0] = 0
    for tsect in range(4):        
        print mon, tsect
        tu, his_pos_trans = pyroms.tools.section_transport_z(u_p, pos_his_v , GRD, xi_rhos[tsect,1], xi_rhos[tsect,0], eta_rhos[tsect], eta_rhos[tsect], h1=0, h2=200)    
        print 'PAST', his_pos_trans

        tu, his_neg_trans = pyroms.tools.section_transport_z(u_p, neg_his_v , GRD, xi_rhos[tsect,1], xi_rhos[tsect,0], eta_rhos[tsect], eta_rhos[tsect], h1=0, h2=200) 
        print 'PAST', his_neg_trans

        tu, fut_pos_trans = pyroms.tools.section_transport_z(u_f, pos_fut_v , GRD, xi_rhos[tsect,1], xi_rhos[tsect,0], eta_rhos[tsect], eta_rhos[tsect], h1=0, h2=200)
        print 'FUTURE', fut_pos_trans

        tu, fut_neg_trans = pyroms.tools.section_transport_z(u_f, neg_fut_v , GRD, xi_rhos[tsect,1], xi_rhos[tsect,0], eta_rhos[tsect], eta_rhos[tsect], h1=0, h2=200)
        print 'FUTURE', fut_neg_trans

        pos_trans_dif[tsect,mon]=(fut_pos_trans-his_pos_trans)/(1000000)
        neg_trans_dif[tsect,mon]=-1*(fut_neg_trans-his_neg_trans)/(1000000)

fig, ax = plt.subplots(4, sharex=True, figsize=(4,6))
fig.subplots_adjust(hspace=.2)
for nt in range(4):
        ax[3-nt].plot([-1,13],[0,0],'-.k')
        if nt == 0:
           ax[3-nt].set_ylim(-3,4)
        else:
           ax[3-nt].set_ylim(-1,3)

        ax[3-nt].set_xlim(-.2,11.2) 
        ax[3-nt].plot(range(12),  pos_trans_dif[nt,:],'-o', mec='k',linewidth=3,markersize=5)
        ax[3-nt].plot(range(12),  neg_trans_dif[nt,:],'-o', mec='k',linewidth=3,markersize=5)
        ax[3-nt].set_xticks([1,4,7,10]) 
        ax[3-nt].set_xticklabels([])
plt.savefig('TRANSPORTS')
plt.show()
