import numpy as np
import netCDF4 as nc
import regrid_atmos 
import pyroms
import matplotlib.pyplot as plt

def save_nc_fil(nf,ncfil,fid,Fout):
    print var[nf]
    fid2 = nc.Dataset(ncfil,'w')

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
    
    fid2.createVariable(var[nf],'f8',('time','eta_rho','xi_rho'),fill_value = np.float(1.0e15))
    fid2.variables['lon'].long_name = fid.variables['lon'].long_name
    fid2.variables[var[nf]].units = fid.variables[var[nf]].units
    fid2.variables[var[nf]].coordinates = fid.variables[var[nf]].coordinates
    fid2.variables[var[nf]].time = fid.variables[var[nf]].time
    fid2.variables[var[nf]][:]=Fout

    fid2.close()

# rid = nc.Dataset('rain_test.nc')
# rain = rid.variables['rain'][:].squeeze()

# CCS details
grd = pyroms.grid.get_ROMS_grid('CCS')

Xout = np.asfortranarray(grd.hgrid.lon_rho.astype(float))
Yout = np.asfortranarray(grd.hgrid.lat_rho.astype(float))

angler = np.asfortranarray(grd.hgrid.angle_rho.astype(float))

Jmax, Imax = Xout.shape

# ERAint files
ncfil =  ['drowned_ERAi_msl_1981-2010_monthly_clim.nc']#,\
#          'drowned_ERAi_precip_1981-2010_monthly_clim.nc',\
#          'drowned_ERAi_q2_1981-2010_monthly_clim.nc',\
#          'drowned_ERAi_radlw_1981-2010_monthly_clim.nc',\
#          'drowned_ERAi_radsw_1981-2010_monthly_clim.nc',\
#          'drowned_ERAi_t2_1981-2010_monthly_clim.nc',\
#          'drowned_ERAi_u10_1981-2010_monthly_clim.nc',\
#          'drowned_ERAi_v10_1981-2010_monthly_clim.nc']

#var = ['Pair','rain','Qair','lwrad_down','swrad','Tair','Uwind','Vwind']
var = ['Pair']

ncdir = '/Users/elizabethdrenkard/external_data/ERAinterim/drowned/80_x_80/'
#ncdir = '/Users/liz.drenkard/external_data/ERAinterim/drowned/80_x_70/'

#for nf in range(len(var)-2):
for nf in range(1):
    file_name = ncdir + ncfil[nf]
    fid = nc.Dataset(file_name)

    Finp = fid.variables[var[nf]][:].squeeze()

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
    Jout = np.zeros((12,Jmax,Imax))
    Iout = np.zeros((12,Jmax,Imax))
    Fout = np.zeros((12,Jmax,Imax))
    for nt in range(Finp.shape[0]):

        Jout[nt,:], Iout[nt,:], Fout[nt,:] = regrid_atmos.regrid_atmos(Xinp, Yinp, np.asfortranarray(Finp[nt,:].squeeze().astype(float)), Amin, Amax, Xout, Yout)
        
        #print np.max(Jout[nt,:].squeeze()), np.min(Jout[nt,:].squeeze())
        #print np.max(Iout[nt,:].squeeze()), np.min(Iout[nt,:].squeeze())
        #plt.pcolor(Jout[nt,:].squeeze())
        #plt.colorbar()
        #plt.show()
    # Save new regridded file
    ncfil2 = 'regridded_' + ncfil[nf]
    save_nc_fil(nf,ncfil2, fid, Fout)

fid.close()

# WINDS A BIT DIFFERENT

#ufil = ncdir + ncfil[-2]
#vfil = ncdir + ncfil[-1]

#ufid = nc.Dataset(ufil)
#vfid = nc.Dataset(vfil)

#FUinp = ufid.variables[var[-2]][:].squeeze()
#FVinp = vfid.variables[var[-1]][:].squeeze()

#AUmin = np.min(FUinp)
#AUmax = np.max(FUinp)
#AVmin = np.max(FVinp)
#AVmax = np.min(FVinp)

#RUN interp routine
#FUout = np.zeros((12,Jmax,Imax))
#FVout = np.zeros((12,Jmax,Imax))

#for nt in range(12):
#    FUout[nt,:], FVout[nt,:] = regrid_atmos.regrid_winds(Xinp, Yinp,    \
#                 np.asfortranarray(FUinp[nt,:].squeeze().astype(float)),\
#                 np.asfortranarray(FVinp[nt,:].squeeze().astype(float)),\
#                 AUmin, AUmax, AVmin, AVmax, angler, Xout, Yout)

# Save new regridded files
#ncfil2 = 'regridded_' + ncfil[-2]
#save_nc_fil(-2,ncfil2, ufid, FUout)  

#ncfil2 = 'regridded_' + ncfil[-1]
#save_nc_fil(-1,ncfil2, vfid, FVout) 

#ufid.close()
#vfid.close()  
