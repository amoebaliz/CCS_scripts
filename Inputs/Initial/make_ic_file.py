import subprocess
import os
import sys
import commands
import numpy as np
import pyroms
import pyroms_toolbox

from remap import remap
from remap_uv import remap_uv

filein = '/glade/p/work/edrenkar/Inputs/construction/Initial/soda3.3.1_5dy_ocean_or_1980_01_03.nc.sub'
dst_dir= './'

tag = '1980_01_03'

print 'Build IC file from the following file:'
print filein
print ' '

# load grids
soda_grd_file = '/glade/p/work/edrenkar/external_data/SODA/SODA3_0.25deg_grid.nc'
src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(soda_grd_file,name='SODA3_0.25deg',xrange=(1,198),yrange=(1,198))
dst_grd = pyroms.grid.get_ROMS_grid('CCS')

zeta_dst_file = dst_dir + dst_grd.name + '_ic_zeta_' + tag + '_' + src_grd.name + '.nc'
temp_dst_file = dst_dir + dst_grd.name + '_ic_temp_' + tag + '_' + src_grd.name + '.nc'
salt_dst_file = dst_dir + dst_grd.name + '_ic_salt_' + tag + '_' + src_grd.name + '.nc'
u_dst_file    = dst_dir + dst_grd.name + '_ic_u_'    + tag + '_' + src_grd.name + '.nc'
v_dst_file    = dst_dir + dst_grd.name + '_ic_v_'    + tag + '_' + src_grd.name + '.nc'

# remap
print filein 
zeta = remap('ssh', filein, src_grd, dst_grd, zeta_dst_file, dst_dir=dst_dir)

# reload grid with zeta (more accurate)
dst_grd = pyroms.grid.get_ROMS_grid('CCS', zeta=zeta)

# regrid temp, salt and velocities
remap('temp',filein, src_grd, dst_grd, temp_dst_file, dst_dir=dst_dir)
remap('salt',filein, src_grd, dst_grd, salt_dst_file, dst_dir=dst_dir)
remap_uv(filein, src_grd, dst_grd, u_dst_file, v_dst_file, dst_dir=dst_dir)

# merge file
ic_file = dst_dir + dst_grd.name + '_ic_' + tag + '_' + src_grd.name + '.nc'

command1 = 'mv '      + zeta_dst_file + ' '    + ic_file
command2 = 'ncks -A ' + temp_dst_file + ' -o ' + ic_file
command3 = 'ncks -A ' + salt_dst_file + ' -o ' + ic_file
command4 = 'ncks -A ' + u_dst_file    + ' -o ' + ic_file
command5 = 'ncks -A ' + v_dst_file    + ' -o ' + ic_file

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

