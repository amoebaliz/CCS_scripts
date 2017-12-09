import subprocess
import os
import sys
import commands
import numpy as np

#increase the maximum number of open files allowed
#import resource
#resource.setrlimit(resource.RLIMIT_NOFILE, (3000,-1))

import pyroms
import pyroms_toolbox

from remap import remap
from remap_uv import remap_uv

#my_year=int(sys.argv[-1])

data_dir = '/glade/p/work/edrenkar/external_data/SODA/'
dst_dir='clim/'

filelst = subprocess.check_output(['ls', data_dir]).replace('/n',' ').split()

#filelstyear = []
#for ff in filelst:
#	if ff.find(str(my_year)) > 0:
#		filelstyear.append(ff)

soda_grd = '/glade/p/work/edrenkar/external_data/SODA/SODA3.4.1_0.25deg_grid_sub.nc'
src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(soda_grd, 'SODA3_0.25deg',xrange=(1,198),yrange=(1,198))
dst_grd = pyroms.grid.get_ROMS_grid('CCS2')

for filein in filelst:
    tag=filein.replace('soda3.4.1_1981-2010_','').replace('.nc','')
    print '\nBuild OBC file for time %s' %filein
    zeta_dst_file = dst_dir + dst_grd.name + '_clim_zeta_' + tag + '_' + src_grd.name + '.nc'
    temp_dst_file = dst_dir + dst_grd.name + '_clim_temp_' + tag + '_' + src_grd.name + '.nc'
    salt_dst_file = dst_dir + dst_grd.name + '_clim_salt_' + tag + '_' + src_grd.name + '.nc'
    u_dst_file    = dst_dir + dst_grd.name + '_clim_u_'    + tag + '_' + src_grd.name + '.nc'
    v_dst_file    = dst_dir + dst_grd.name + '_clim_v_'    + tag + '_' + src_grd.name + '.nc'

    # remap ssh
    zeta = remap('ssh', data_dir + filein, src_grd, dst_grd, zeta_dst_file, dst_dir=dst_dir)

    # reload grid with zeta (more accurate)
    dst_grd = pyroms.grid.get_ROMS_grid('CCS2', zeta=zeta)

    # regrid temp, salt and velocities
    remap('temp',data_dir + filein, src_grd, dst_grd, temp_dst_file, dst_dir=dst_dir)
    remap('salt',data_dir + filein, src_grd, dst_grd, salt_dst_file, dst_dir=dst_dir)
    remap_uv(data_dir + filein, src_grd, dst_grd, u_dst_file, v_dst_file, dst_dir=dst_dir)

    # merge file
    clim_file = dst_dir + dst_grd.name + '_clim_' + tag + '_' + src_grd.name + '.nc'

    command1 = 'mv '      + zeta_dst_file + ' '    + clim_file
    command2 = 'ncks -A ' + temp_dst_file + ' -o ' + clim_file
    command3 = 'ncks -A ' + salt_dst_file + ' -o ' + clim_file
    command4 = 'ncks -A ' + u_dst_file    + ' -o ' + clim_file
    command5 = 'ncks -A ' + v_dst_file    + ' -o ' + clim_file

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
