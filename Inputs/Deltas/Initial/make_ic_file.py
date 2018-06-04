import subprocess
import os
import sys
import commands
import numpy as np
import pyroms
import pyroms_toolbox

from remap import remap
from remap_uv import remap_uv

#start = datetime.now()
dst_dir='init/'
for run in [17]:
    runstr = str(run).zfill(3)
    filein = '/Users/elizabethdrenkard/Desktop/LENS_deltas_ocn/' + runstr + '/' + runstr + '_ocn_01.nc' 
    print 'Build IC file from the following file:'
    print filein
    print ' '

    # load grids
    soda_grd = '/glade/p/work/edrenkar/external_data/LENS/LENS_grid.nc'
    soda_grd = '/Users/elizabethdrenkard/Desktop/LENS_deltas_ocn/017/LENS_grid.nc'
    src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(soda_grd, 'LENS',xrange=(1,78),yrange=(1,78))
    dst_grd = pyroms.grid.get_ROMS_grid('CCS')

    tag = runstr

    zeta_dst_file = dst_dir + dst_grd.name + '_ic_zeta_' + tag + '_' + src_grd.name + '.nc'
    temp_dst_file = dst_dir + dst_grd.name + '_ic_temp_' + tag + '_' + src_grd.name + '.nc'
    salt_dst_file = dst_dir + dst_grd.name + '_ic_salt_' + tag + '_' + src_grd.name + '.nc'
    u_dst_file    = dst_dir + dst_grd.name + '_ic_u_'    + tag + '_' + src_grd.name + '.nc'
    v_dst_file    = dst_dir + dst_grd.name + '_ic_v_'    + tag + '_' + src_grd.name + '.nc'

    # remap
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

