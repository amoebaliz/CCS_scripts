import numpy as np
import netCDF4 as nc
from os import listdir
from os.path import isfile, join
import glob
import pyroms
import matplotlib.pyplot as plt


def make_ncfil(var,var_str):

    ncfil = '/Users/liz.drenkard/Desktop/monthly_CCS_slices/' + \
            var_str + '_' + str(depth) + 'm_his_' + str(mn+1).zfill(2) + '_' + str(yr).zfill(4) + '.nc'
    
    if var_str in ['temp', 'salt']:
       Cpos = 'rho'
    else:
       Cpos = var_str

    eta_dim = 'eta_' + Cpos
    xi_dim = 'xi_'+ Cpos
    eta_str = 'lat_' + Cpos + '.shape[0]'
    xi_str = 'lat_' + Cpos + '.shape[1]'
    
    fid2 = nc.Dataset(ncfil,'w')
    fid2.createDimension('ocean_time', None)    
    fid2.createDimension(eta_dim, eval(eta_str))
    fid2.createDimension(xi_dim, eval(xi_str))


    fid2.createVariable('ocean_time', 'f8', ('ocean_time'))
    fid2.variables['ocean_time'] = fid.variables['ocean_time'] 
    
    print fid.variables[var_str]
    fid2.createVariable(var_str,'f8',('ocean_time',eta_dim, xi_dim))  
    if var_str != 'salt':
       fid2.variables[var_str].units = fid.variables[var_str].units
    fid2.variables[var_str].long_name = fid.variables[var_str].long_name
    fid2.variables[var_str][0,:]= var

    fid2.close()

# ROMS Grid information
grd = pyroms.grid.get_ROMS_grid('CCS')
lat_rho = grd.hgrid.lat_rho[:]
lat_u = grd.hgrid.lat_u[:]
lat_v = grd.hgrid.lat_v[:]

# Depth of interest meters
depth=100

for yr in range(3,12+1):
    #print yr

    for mn in range(12):
        #print mn+1
        mypath = '/Volumes/Abalone/CCS/his2/' + str(yr).zfill(4) + \
                 '/CCS-LD.HCo02Y_avg_' + str(yr).zfill(4) + \
                 '-' + str(mn+1).zfill(2) + '*T12:00:00.nc'
           
        # list all files of given year-month make sure they exist. 
        # NOTE: day order in month will not be preserved BUT, this 
        #       should be unnecessary if commputing monthly means
        n=0
        for f in glob.glob(mypath):
            print f 
            fid = nc.Dataset(f)
            #zeta = fid.variables['zeta'][:].squeeze() 
            temp = fid.variables['temp'][:].squeeze()
            salt = fid.variables['salt'][:].squeeze()
            u = fid.variables['u'][:].squeeze()
            v = fid.variables['v'][:].squeeze()
            
            # Initalize storage fields
            if n==0:
               t_stor = np.ma.zeros(lat_rho.shape)
               s_stor = np.ma.zeros(lat_rho.shape)  
               u_stor = np.ma.zeros(lat_u.shape)
               v_stor = np.ma.zeros(lat_v.shape)

            # UPDATE STORAGE FIELDS WITH CALCULATED VALUE AT DESIRED DEPTH
            # NOTE: Using default linear interpolation (as opposed to a spine)
            t_stor+= pyroms.tools.zslice(temp, depth, grd,  Cpos='rho')[0]
            s_stor+= pyroms.tools.zslice(salt, depth, grd,  Cpos='rho')[0]
            u_stor+= pyroms.tools.zslice(u, depth, grd, Cpos='u')[0]
            v_stor+= pyroms.tools.zslice(v, depth, grd, Cpos='v')[0]

            n+=1

        # Generate monthly mean fields by dividing by n
        t_depth = t_stor/n
        s_depth = s_stor/n
        u_depth = u_stor/n
        v_depth = v_stor/n

        # Save as monthly netCDF file
        make_ncfil(t_depth,'temp')
        make_ncfil(s_depth,'salt')
        make_ncfil(u_depth,'u')
        make_ncfil(v_depth,'v') 

        print str(yr).zfill(4), str(mn+1).zfill(2)
