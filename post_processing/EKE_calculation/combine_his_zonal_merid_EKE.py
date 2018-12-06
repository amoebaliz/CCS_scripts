import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# ROMS Grid information
grdfile = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid = nc.Dataset(grdfile)
# Subset of grid because velocity fields are not averaged to outter rho points
mask = fid.variables['mask_rho'][1:-1,1:-1]
lat = fid.variables['lat_rho'][1:-1,1:-1]
lon = fid.variables['lon_rho'][1:-1,1:-1]

ny, nx = lat.shape

EKE_stor=np.zeros(lat.shape)
for nt in range(10):
    ncfil = 'CCS_MAM_uveke_' +str(nt+3).zfill(4) + '.nc'
    #fid = nc.Dataset(ncfil,'a')
    fid = nc.Dataset(ncfil)  
    u_eke = fid.variables['u_eke'][:].squeeze()
    v_eke = fid.variables['v_eke'][:].squeeze()
    temp  = fid.variables['temp'][:].squeeze()

    # average u over xi dimension
    EKE = ((u_eke[1:-1,:-1]+u_eke[1:-1,1:])/2. + (v_eke[:-1,1:-1]+v_eke[1:,1:-1])/2.)/2
    # Save SST and EKE as netCDF  
    # Create new file    
    new_ncfil =  'CCS_MAM_his_EKE_' + str(nt+3).zfill(4) + '.nc' 
    fid2 = nc.Dataset(new_ncfil, 'w') 
    
    # Write EKE to file
    fid2.createDimension('ocean_time', None)
    fid2.createDimension('eta_rho', ny) 
    fid2.createDimension('xi_rho', nx)  
    
    time = fid2.createVariable('ocean_time', 'f8', ('ocean_time'))
    fid2.variables['ocean_time'].units = 'year_number'
    fid2.variables['ocean_time'][0] = nt+1

    fid2.createVariable('lat_rho','f8',('eta_rho','xi_rho'))
    fid2.variables['lat_rho'].units = 'latitude'
    fid2.variables['lat_rho'][:]=lat

    fid2.createVariable('lon_rho','f8',('eta_rho','xi_rho'))
    fid2.variables['lon_rho'].units = 'longitude'
    fid2.variables['lon_rho'][:]=lon

    SST = fid2.createVariable('SST','f8',('ocean_time','eta_rho','xi_rho'))
    fid2.variables['SST'].units = 'oC'
    fid2.variables['SST']._FillValue=np.max(temp)
    fid2.variables['SST'].long_name = 'Sea Surface Temperature'
    fid2.variables['SST'][0,:]= temp[1:-1,1:-1]
     
    EKE2 = fid2.createVariable('EKE','f8',('ocean_time','eta_rho','xi_rho'))
    fid2.variables['EKE'].units = 'm^2/s^2'
    fid2.variables['EKE']._FillValue=np.max(EKE)
    fid2.variables['EKE'].long_name = 'total eddy kinetic energy'
    fid2.variables['EKE'][0,:]= EKE

    fid2.close()


    EKE_stor = EKE_stor + EKE
    print np.sum(EKE) 
    #plt.figure()
    #plt.pcolormesh(lon,lat,EKE,vmin=0,vmax=0.08)
    #plt.colorbar()
    plt.plot(nt,np.sum(EKE),'-o')
    fid.close()     

plt.figure()
plt.pcolormesh(lon,lat,EKE_stor/10.,vmin=0,vmax=0.05)
plt.colorbar()
#plt.show()



