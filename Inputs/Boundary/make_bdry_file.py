import subprocess
import os
import sys
import commands
import numpy as np
import pyroms
import pyroms_toolbox

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv

soda_grd = '/glade/p/work/edrenkar/external_data/SODA/SODA3.4.1_0.25deg_grid.nc'

src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(soda_grd, name='SODA3_0.25deg',xrange=(1,198),yrange=(1,198))
dst_grd = pyroms.grid.get_ROMS_grid('CCS')

data_dir = '/glade/p/work/edrenkar/external_data/SODA/'
dst_dir='bdry/'

lst_file = []
lst = []

for run in range(12):
    lst.append(data_dir + 'soda3.4.1_1981-2010_clim_' + str(run+1).zfill(2) + '.nc') 

lst_file = lst_file + lst

for filein in lst_file:
        tag=filein.replace(data_dir + 'soda3.4.1_1981-2010_clim_','').replace('.nc','')
        print '\nBuild OBC file for time %s' %filein
        zeta_dst_file = dst_dir + dst_grd.name + '_bdry_zeta_' + tag + '_' + src_grd.name + '.nc'
        temp_dst_file = dst_dir + dst_grd.name + '_bdry_temp_' + tag + '_' + src_grd.name + '.nc'
        salt_dst_file = dst_dir + dst_grd.name + '_bdry_salt_' + tag + '_' + src_grd.name + '.nc'
        u_dst_file    = dst_dir + dst_grd.name + '_bdry_u_'    + tag + '_' + src_grd.name + '.nc'
        v_dst_file    = dst_dir + dst_grd.name + '_bdry_v_'    + tag + '_' + src_grd.name + '.nc'

        # remap ssh
        zeta = remap_bdry('ssh', filein, src_grd, dst_grd, zeta_dst_file, dst_dir=dst_dir)

        # reload grid with zeta (more accurate)
        dst_grd = pyroms.grid.get_ROMS_grid('CCS', zeta=zeta)

        # regrid temp, salt and velocities
        remap_bdry('temp', filein, src_grd, dst_grd, temp_dst_file, dst_dir=dst_dir)
        remap_bdry('salt', filein, src_grd, dst_grd, salt_dst_file, dst_dir=dst_dir)
        remap_bdry_uv(filein, src_grd, dst_grd, u_dst_file, v_dst_file, dst_dir=dst_dir)

        # merge file
        bdry_file = dst_dir + dst_grd.name + '_bdry_' + tag + '_' + src_grd.name + '.nc'

        command1 = 'mv '      + zeta_dst_file + ' '    + bdry_file
        command2 = 'ncks -A ' + temp_dst_file + ' -o ' + bdry_file
        command3 = 'ncks -A ' + salt_dst_file + ' -o ' + bdry_file
        command4 = 'ncks -A ' + u_dst_file    + ' -o ' + bdry_file
        command5 = 'ncks -A ' + v_dst_file    + ' -o ' + bdry_file

        subprocess.call(command1,shell=True)
        subprocess.call(command2,shell=True)
        subprocess.call(command3,shell=True)
        subprocess.call(command4,shell=True)
        subprocess.call(command5,shell=True)

        # clean up
        os.remove(temp_dst_file)
        os.remove(salt_dst_file)
        os.remove(u_dst_file)
        os.remove(v_dst_file)
