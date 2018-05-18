import subprocess
import os
import sys 

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


dfs_dir = '/Users/liz.drenkard/external_data/DFS5.2/climatologies/'
era_dir = '/Users/liz.drenkard/external_data/ERAinterim/Forecast/climatologies/'

#filelst_dfs = subprocess.check_output(['ls', dfs_dir]).replace('/n',' ').split()
#filelst_era = subprocess.check_output(['ls', era_dir]).replace('/n',' ').split()

filelst = subprocess.check_output(['ls', dfs_dir]).replace('/n',' ').split()

varlst = ['Pair', 'rain', 'Qair','lwrad_down','swrad','Tair','Uwind', 'Vwind']
# CHECK SAME FILE NAMES
#for (d_fil,e_fil) in  zip(filelst_dfs,filelst_era): 
#    if d_fil == e_fil:
#       print 'MEEEEEP'

for (fil,var) in zip(filelst,varlst):
    print var
    ncfil_d = dfs_dir + fil
    ncfil_e = era_dir + fil

    fid_d = nc.Dataset(ncfil_d)
    fid_e = nc.Dataset(ncfil_e)

    var_d = fid_d.variables[var][:]
    var_e = fid_e.variables[var][:]
    
    for nt in range(var_d.shape[0]):
        plt.figure()
        plt.pcolor(var_d[nt,:].squeeze()-var_e[nt,:].squeeze(),vmin=-10,vmax=10)
        plt.colorbar()

    plt.show()





#plt.show()




