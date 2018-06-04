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

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv

#my_year=int(sys.argv[-1])

dst_dir='bdry/'
for run in [17]:
    runstr = str(run).zfill(3)
    data_dir = '/Users/elizabethdrenkard/Desktop/LENS_deltas_ocn/' + runstr + '/'
    filelst = subprocess.check_output(['ls', data_dir]).replace('/n',' ').split()
    print filelst
    soda_grd = '/glade/p/work/edrenkar/external_data/LENS/LENS_grid.nc'
    soda_grd = '/Users/elizabethdrenkard/Desktop/LENS_deltas_ocn/017/LENS_grid.nc'
    src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(soda_grd, 'LENS',xrange=(1,78),yrange=(1,78))
    dst_grd = pyroms.grid.get_ROMS_grid('CCS')

    for filein in filelst:
        tag=filein.replace('LENS_DELTA_','').replace('.nc','')
        print '\nBuild OBC file for time %s' %filein
        print data_dir + filein
        zeta_dst_file = dst_dir + dst_grd.name + '_bdry_zeta_' + tag + '_' + src_grd.name + '.nc'
        temp_dst_file = dst_dir + dst_grd.name + '_bdry_temp_' + tag + '_' + src_grd.name + '.nc'
        salt_dst_file = dst_dir + dst_grd.name + '_bdry_salt_' + tag + '_' + src_grd.name + '.nc'
        u_dst_file    = dst_dir + dst_grd.name + '_bdry_u_'    + tag + '_' + src_grd.name + '.nc'
        v_dst_file    = dst_dir + dst_grd.name + '_bdry_v_'    + tag + '_' + src_grd.name + '.nc'

        # remap ssh
        zeta = remap_bdry('ssh', data_dir + filein, src_grd, dst_grd, zeta_dst_file, dst_dir=dst_dir)

        # reload grid with zeta (more accurate)
        dst_grd = pyroms.grid.get_ROMS_grid('CCS', zeta=zeta)

        # regrid temp, salt and velocities
        remap_bdry('temp',data_dir + filein, src_grd, dst_grd, temp_dst_file, dst_dir=dst_dir)
        remap_bdry('salt',data_dir + filein, src_grd, dst_grd, salt_dst_file, dst_dir=dst_dir)
        remap_bdry_uv(data_dir + filein, src_grd, dst_grd, u_dst_file, v_dst_file, dst_dir=dst_dir)

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
