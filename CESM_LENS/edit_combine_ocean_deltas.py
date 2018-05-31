import numpy as np
import netCDF4 as nc

def build_ocean_fil(run,mon):
    runstr = str(run).zfill(3)
    monstr = str(mon).zfill(2)
          
    ncfil = '/glade/p/work/edrenkar/Inputs/construct/Deltas/LENS_deltas_ocn/' + runstr + '/' + runstr + '_ocn_' + monstr + '.nc'

    fid = nc.Dataset(ncfil,'w')

    # Dimensions 
    fid.createDimension('time', None)
    fid.createDimension('st_ocean',60)
    fid.createDimension('sw_ocean',60)
    fid.createDimension('yt_ocean',80)
    fid.createDimension('yu_ocean',80)
    fid.createDimension('xt_ocean',80)
    fid.createDimension('xu_ocean',80)

    # TIME
    src_fil = dir + runstr + '_SSH_clim_delta.nc'
    fid.createVariable('time', 'f8', ('time'))
    fid.variables['time'].units = "days since 1900-01-01 00:00:00"
    fid.variables['time'][:] = nc.Dataset(src_fil).variables['time'][mon-1].squeeze()

    # SSH
    src_fil = dir + runstr + '_SSH_clim_delta.nc'
    srcSSH = nc.Dataset(src_fil).variables['SSH']
    fid.createVariable('ssh', 'f8', ('time','yt_ocean', 'xt_ocean'),fill_value = np.float(1.0e15))
    fid.variables['ssh'].missing_value = np.float(1.0e15)
    fid.variables['ssh'].units = "meters"
    fid.variables['ssh'][0] = srcSSH[mon-1,:].squeeze()/100.

    # TEMP
    src_fil = dir + runstr + '_TEMP_clim_delta.nc'
    srcTEMP = nc.Dataset(src_fil).variables['TEMP']
    fid.createVariable('temp', 'f8', ('time','st_ocean','yt_ocean', 'xt_ocean'),fill_value = np.float(1.0e15))
    fid.variables['temp'].missing_value = np.float(1.0e15)
    fid.variables['temp'].units = srcTEMP.units
    fid.variables['temp'][0,:] = srcTEMP[mon-1,:].squeeze()

    # SALT
    src_fil = dir + runstr + '_SALT_clim_delta.nc'
    srcSALT = nc.Dataset(src_fil).variables['SALT']
    fid.createVariable('salt', 'f8', ('time','st_ocean','yt_ocean', 'xt_ocean'),fill_value = np.float(1.0e15))
    fid.variables['salt'].missing_value = np.float(1.0e15)
    fid.variables['salt'].units = srcSALT.units
    fid.variables['salt'][0,:] = srcSALT[mon-1,:].squeeze()

    # UVEL
    src_fil = dir + runstr + '_UVEL_clim_delta.nc'
    srcUVEL = nc.Dataset(src_fil).variables['UVEL']
    fid.createVariable('u', 'f8', ('time','st_ocean','yu_ocean', 'xu_ocean'),fill_value = np.float(1.0e15))
    fid.variables['u'].missing_value = np.float(1.0e15)
    fid.variables['u'].units = "meters/s"
    fid.variables['u'][0,:] = srcUVEL[mon-1,:].squeeze()/100.

    # VVEL
    src_fil = dir + runstr + '_VVEL_clim_delta.nc'
    srcVVEL = nc.Dataset(src_fil).variables['VVEL']
    fid.createVariable('v', 'f8', ('time','st_ocean','yu_ocean', 'xu_ocean'),fill_value = np.float(1.0e15))
    fid.variables['v'].missing_value = np.float(1.0e15)
    fid.variables['v'].units = "meters/s"
    fid.variables['v'][0,:] = srcVVEL[mon-1,:].squeeze()/100.

    fid.close() 
    print ncfil
# -----------
vars={'SSH','TEMP','SALT','UVEL','VVEL'}
new_vars = {'ssh', 'temp','salt','u','v'}

dir = '/glade/p/work/edrenkar/external_data/LENS/difs/'
dir = '/glade/u/home/edrenkar/TOOLS/CCS_scripts/CESM_LENS/difs/'
grd = nc.Dataset('/glade/p/work/edrenkar/external_data/LENS/LENS_grid.nc')

for run in [17]:

    for mon in range(12):

        # CREATE NEW MONTHLY CLIM FILE WITH ALL OCEAN FIELDS
        build_ocean_fil(run,mon+1)



