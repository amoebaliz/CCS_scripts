import matplotlib
matplotlib.use('Agg')
import numpy as np
import netCDF4 as nc
from datetime import datetime
import sys

invar = 'salt'
outvar = 'SSS'
outtime = 'sss_time'

#read grid and variable attributes from the first file
spval = np.float(-99)
units = []
long_name = []

ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
dtot=0
newtime=np.zeros(len(ndays))

# establish sss_time values for cycling over 365d yr
# (i.e., middle of Jan, Feb, Mar ... Dec)
for nmon in range(len(ndays)):
    newtime[nmon] = dtot + ndays[nmon]/2.0
    dtot+=ndays[nmon]

# read in original salinity file
ncfil = 'phc3.0_monthly.nc'
fid = nc.Dataset(ncfil)
lon = fid.variables['lon'][:]
lat = fid.variables['lat'][:]

# create new SSS file with information from the old file 
outfile = 'sss_monthly_climatology.nc'
newnc = nc.Dataset(outfile, 'w', format='NETCDF3_CLASSIC')
newnc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
newnc.conventions = 'CF-1.6'
newnc.title = 'A Global Ocean Hydrography with a High Quality Arctic Ocean Version 3.0'
newnc.reference = 'Steele, Morley and Ermold (2001), PHC: A global ocean hydrography with a high-quality Arctic Ocean, J. Clim'
newnc.source = 'Polar Science Center/Applied Physics Lab/University of Washington'
newnc.institution = 'University of Washington'

newnc.createDimension('lon', np.size(lon))
newnc.createDimension('lat', np.size(lat))
newnc.createDimension(outtime, None)

newnc.createVariable('lon', 'f8', ('lon'))
newnc.variables['lon'].long_name = 'longitude'
newnc.variables['lon'].units = 'degrees_east'
newnc.variables['lon'][:] = lon

newnc.createVariable('lat', 'f8', ('lat'))
newnc.variables['lat'].long_name = 'latitude'
newnc.variables['lat'].units = 'degrees_north'
newnc.variables['lat'][:] = lat

newnc.createVariable(outtime, 'f8', (outtime))
newnc.variables[outtime].units = 'Days'
newnc.variables[outtime][:] = newtime

newnc.createVariable(outvar, 'f', (outtime, 'lat', 'lon'), fill_value=spval)
newnc.variables[outvar].missing_value = spval
newnc.variables[outvar].long_name = 'Salinity, modified 10m average'
newnc.variables[outvar].units = 'psu'
newnc.variables[outvar].coordinates = 'lon lat'

# get data
var = np.mean(fid.variables[invar][:,0:2,:,:],axis=1).squeeze()

# for PCH3.0 method
fid2 = nc.Dataset('sss_monthly_climatology_flooded_first.nc')
salt_0 = fid2.variables['SSS'][:]

ztmp = np.empty(var.shape)
### flood land values ###
# following algorithms from Raf's creeping sea fortran algorithm #
nmon,ny,nx = var.shape
for jmon in range(nmon):
    a = var[jmon,:,:].copy()
    salt_sub = salt_0[jmon,:,:].copy()
    # land sea weightings
    sor   = np.zeros((ny,nx))
    sor[np.isnan(a)]=.6

    # guess for land values - avg salt value
    # a[np.isnan(a)] = 34
    a[np.isnan(a)] = salt_sub[np.isnan(a)]   
 
    nt = 0
    nmax=40
    resmax=0.0
    ct=1.e-5
    absres=0.0
    while (nt < nmax) and (absres < ct):
#    while (nt < 2):
      res_val = np.zeros((ny,nx))
      # switched dimensional naming convention
      for i in range(1,ny-1):           
          for j in range(1,nx-1):
              res_val[i,j] = .25*(a[i-1,j] + a[i+1,j] + a[i,j-1] + a[i,j+1])-a[i,j]
      for i in range(1,ny-1):
          for j in range(1,nx-1):
              res_val[i,j] = res_val[i,j]*sor[i,j] 
              a[i,j] = a[i,j] + res_val[i,j]
              
              absres = abs(res_val[i,j])
              resmax = max(absres,resmax)
 
      #NO FLUX EDGE VALUES
      #E/W edge vales
      for i in range(ny):
          a[i,0]  = a[i,1]
          a[i,nx-1] = a[i,nx-2]
      #N/S edge values 
      for j in range(nx):
          a[0,j]  = a[1,j]
          a[ny-1,j] = a[ny-2,j]
           
      nt += 1
      # update mask AFTER for looks above so as not to bias the 
      # the flood values to the for loop origin
      #zmask[zmask_sub==1]=1
      #zmask_halo2pt[zmask_halo2pt_sub==1]=1

    ztmp[jmon,:,:] = a
newnc.variables[outvar][:] = ztmp
newnc.close()
fid.close()

