import numpy as np
import netCDF4 as nc
import regrid_atmos 
import pyroms
import matplotlib.pyplot as plt

# rid = nc.Dataset('rain_test.nc')
# rain = rid.variables['rain'][:].squeeze()

# CCS details
grd = pyroms.grid.get_ROMS_grid('CCS')

Xout = np.asfortranarray(grd.hgrid.lon_rho.astype(float))
Yout = np.asfortranarray(grd.hgrid.lat_rho.astype(float))

Jmax, Imax = Xout.shape

# ERAint files
ncfil =  ['drowned_ERAi_msl_1981-2010_monthly_clim.nc',\
          'drowned_ERAi_precip_1981-2010_monthly_clim.nc',\
          'drowned_ERAi_q2_1981-2010_monthly_clim.nc',\
          'drowned_ERAi_radlw_1981-2010_monthly_clim.nc',\
          'drowned_ERAi_radsw_1981-2010_monthly_clim.nc',\
          'drowned_ERAi_t2_1981-2010_monthly_clim.nc',\
          'drowned_ERAi_u10_1981-2010_monthly_clim.nc',\
          'drowned_ERAi_v10_1981-2010_monthly_clim.nc']

#ncfil = ['drowned_ERAi_radsw_1981-2010_monthly_clim.nc']

var = ['Pair','rain','Qair','lwrad_down','swrad','Tair','Uwind','Vwind']

#var = ['swrad']
for nfil in range(len(ncfil)):

    file_name= '/Users/elizabethdrenkard/external_data/ERAinterim/drowned/' + ncfil[nfil]

    #fid = nc.Dataset(ncfil[nfil])
    fid = nc.Dataset(file_name)

    Finp = fid.variables[var[nfil]][:].squeeze()

    lat = fid.variables['lat'][:].squeeze().astype(float)
    lon = fid.variables['lon'][:].squeeze().astype(float)
    lon[lon>180]=lon[lon>180]-360

    Ny = len(lat)
    Nx = len(lon)
    
    Yinp = np.asfortranarray(np.transpose(np.tile(lat,(Nx,1))))
    Xinp = np.asfortranarray(np.tile(lon,(Ny,1)))

    Amin = np.min(Finp)
    Amax = np.max(Finp)

    # Run interp routine
    Fout = np.zeros((12,Jmax,Imax))
    for nt in range(Finp.shape[0]):

        Fout[nt,:] = regrid_atmos.regrid_atmos(np.asfortranarray(Finp[nt,:].squeeze().astype(float)), Xinp, Yinp, Amin, Amax, Xout, Yout)

    # Save new regridded file
    ncfil2 = 'regridded_' + ncfil[nfil]
    fid2 = nc.Dataset(ncfil2,'w')

    fid2.createDimension('time', None)
    fid2.createDimension('eta_rho', Jmax)
    fid2.createDimension('xi_rho', Imax)

    fid2.createVariable('time', 'f8', ('time'))
    fid2.variables['time'].units = fid.variables['time'].units
    fid2.variables['time'].cycle_length = fid.variables['time'].cycle_length
    fid2.variables['time'][:] = fid.variables['time'][:]

    fid2.createVariable('lat','f8',('eta_rho','xi_rho'))
    fid2.variables['lat'].long_name = fid.variables['lat'].long_name
    fid2.variables['lat'].units = fid.variables['lat'].units
    fid2.variables['lat'][:]=Yout

    fid2.createVariable('lon','f8',('eta_rho','xi_rho'))
    fid2.variables['lon'].long_name = fid.variables['lon'].long_name
    fid2.variables['lon'].units = fid.variables['lon'].units
    fid2.variables['lon'][:]=Xout
    
    fid2.createVariable(var[nfil],'f8',('time','eta_rho','xi_rho'),fill_value = np.float(1.0e15))
    fid2.variables['lon'].long_name = fid.variables['lon'].long_name
    fid2.variables[var[nfil]].units = fid.variables[var[nfil]].units
    fid2.variables[var[nfil]].coordinates = fid.variables[var[nfil]].coordinates
    fid2.variables[var[nfil]].time = fid.variables[var[nfil]].time
    fid2.variables[var[nfil]][:]=Fout

    fid2.close()
fid.close()
