import matplotlib
matplotlib.use('Agg')
import numpy as np
import netCDF4 as nc
from datetime import datetime
import sys

def gauss_smooth(a,mask):

    nt = 0
#    nmax=10
    resmax=0.0
    ct=1.e-5
    absres=0.0
    sor   = np.zeros((ny,nx))
    sor[mask==0]=.6

    while (nt < nmax) and (absres < ct):    
          res_val = np.zeros((ny,nx))
          for j in range(1,ny-1):
              for i in range(1,nx-1):
                  res_val[j,i] = .25*(a[j-1,i] + a[j+1,i] + a[j,i-1] + a[j,i+1])-a[j,i]

          for j in range(1,ny-1):
              for i in range(1,nx-1):
                  res_val[j,i] = res_val[j,i]*sor[j,i]
                  a[j,i] = a[j,i] + res_val[j,i]

                  absres = abs(res_val[j,i])
                  resmax = max(absres,resmax)

          #NO FLUX EDGE VALUES
          #E/W edge vales
          for j in range(ny):
              a[j,0]  = a[j,1]
              a[j,nx-1] = a[j,nx-2]
          #N/S edge values 
          for i in range(nx):
              a[0,i]  = a[1,j]
              a[ny-1,i] = a[ny-2,i]

          nt += 1 
    return a


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
ztmp = np.empty(var.shape)
### flood land values ###
# following algorithms from Raf's creeping sea fortran algorithm #
nmon,ny,nx = var.shape
for jmon in range(nmon):
    nedit=0
    var_tmp = var[jmon,:,:].copy()

    zmask   = np.ones((ny,nx))
    zmask[np.isnan(var_tmp)]=0
    orig_mask = zmask.copy()    

    zmask_halo2pt = np.zeros((ny+2,nx+2))
    zmask_halo2pt[1:ny+1,1:nx+1] = zmask

    ztmp_halo2pt  = np.ones((ny+2,nx+2))*spval
    ztmp_halo2pt[1:ny+1,1:nx+1]  = var_tmp
    ztmp_halo2pt[np.isnan(ztmp_halo2pt)]=spval

    # Halo arrangement around central point
    # for which generating weighted mean sea value:

    # c1 -- c2 -- c3
    # |           |
    # c4          c5
    # |           |
    # c6 -- c7 -- c8

    nt = 0
    nmax=500
    while np.any(zmask==0) and (nt < nmax):
      zmask_sub=np.zeros((ny,nx))
      zmask_halo2pt_sub=np.zeros((ny+2,nx+2))
      for jj in range(ny):
          for ji in range(nx):
              # compute indexes surrounding the point at [jj,ji]
              jjm1 = jj-1 ; jjp1 = jj+1
              jim1 = ji-1 ; jip1 = ji+1
         
              # move up and over to account for halo in storage variables
              jjm1h = jjm1 + 1 ; jjp1h = jjp1 + 1
              jim1h = jim1 + 1 ; jip1h = jip1 + 1
              jih   = ji   + 1 ; jjh   = jj   + 1

              # If dealing with a land point at [jj,ji], we need to fill it
              if (zmask[jj,ji] == 0 ):
                 c6 = 1 * zmask_halo2pt[jjm1h,jim1h]
                 c7 = 2 * zmask_halo2pt[jjm1h,jih]
                 c8 = 1 * zmask_halo2pt[jjm1h,jip1h]

                 c4 = 2 * zmask_halo2pt[jjh,jim1h]
                 c5 = 2 * zmask_halo2pt[jjh,jip1h]

                 c1 = 1 * zmask_halo2pt[jjp1h,jim1h]
                 c2 = 2 * zmask_halo2pt[jjp1h,  jih]
                 c3 = 1 * zmask_halo2pt[jjp1h,jip1h]

                 ctot = c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8

                 if (ctot >= 3 ):
                    nedit+=1
                    # compute the new value for this point
                    zval  = (c6 * ztmp_halo2pt[jjm1h,jim1h] + \
                             c7 * ztmp_halo2pt[jjm1h,  jih] + \
                             c8 * ztmp_halo2pt[jjm1h,jip1h] + \
                             c4 * ztmp_halo2pt[jjh,  jim1h] + \
                             c5 * ztmp_halo2pt[jjh,  jip1h] + \
                             c1 * ztmp_halo2pt[jjp1h,jim1h] + \
                             c2 * ztmp_halo2pt[jjp1h,  jih] + \
                             c3 * ztmp_halo2pt[jjp1h,jip1h])/ \
                             ( ctot )
                    # update value in field array
                    ztmp_halo2pt[jjh,jih] = zval
   
                    # set the mask to sea
                    zmask_sub[jj,ji] = 1
                    zmask_halo2pt_sub[jjh,jih] = 1
      nt += 1
      # update mask AFTER for looks above so as not to bias the 
      # the flood values to the for loop origin
      zmask[zmask_sub==1]=1
      zmask_halo2pt[zmask_halo2pt_sub==1]=1

    ztmp[jmon,:,:] = ztmp_halo2pt[1:ny+1,1:nx+1]
    #ztmp[jmon,:,:] = gauss_smooth(ztmp_halo2pt[1:ny+1,1:nx+1],orig_mask)
    
# ztmp gives us the "best land estimates" for using with pch3.0 gaussian ... thing

newnc.variables[outvar][:] = ztmp
newnc.close()
fid.close()

