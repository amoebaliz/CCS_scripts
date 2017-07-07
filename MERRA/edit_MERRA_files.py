import matplotlib
matplotlib.use('Agg')
import numpy as np
import netCDF4 as nc
from datetime import datetime
import sys

#invar = ['QV2M','SLP','T2M','U2M','V2M','SWGDN','LWGAB','PRECTOT']
invar = ['U10M','V10M']
#outvar = ['Qair','Pair','Tair','Uwind','Vwind','swrad','lwrad_down','rain']
outvar = ['Uwind','Vwind']
fillab = ['slv','slv','slv','slv','slv','rad','rad','flx']
outtime = 'time'

#read grid and variable attributes from the first file
spval = np.empty(len(invar))
units = []
long_name = []

ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
dtot=np.zeros(len(invar))

#for nvar in range(len(invar)):
for nvar in range(2):
    if nvar<5:
       fillab = 'slv_10m'
       #fillab = 'slv'
    elif nvar>6:
       fillab = 'flx'
    else:
       fillab = 'rad'

    for nmon in range(12):
        #ncfil = './drowned/drowned_MERRA2_' + fillab + '_1981-2010_' + str(nmon+1).zfill(2) + '.nc'
        ncfil = 'MERRA2_' + fillab + '_1981-2010_' + str(nmon+1).zfill(2) + '.nc'
        #print ncfil
        fid = nc.Dataset(ncfil)
        if nmon == 0:
           lon = fid.variables['lon'][:]
           lat = fid.variables['lat'][:]
           spval[nvar] = fid.variables[invar[nvar]].missing_value
           if invar[nvar]=='T2M':
              units.append('Celsius')
           else:
              units.append(fid.variables[invar[nvar]].units) 
           long_name.append(fid.variables[invar[nvar]].long_name)
           
        #create ROMS forcing file
        outfile = 'MERRA_' + outvar[nvar] + '_1981-2010_MONTHLY_CLIM_' + str(nmon+1).zfill(2) + '.nc'
        #outfile = 'MERRA_' + outvar[nvar] + '_10M_1981-2010_MONTHLY_CLIM_' + str(nmon+1).zfill(2) + '.nc'
        newnc = nc.Dataset(outfile, 'w', format='NETCDF3_CLASSIC')
        newnc.Author = sys._getframe().f_code.co_name
        newnc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        newnc.title = 'MERRA dataset. Modern Era Retrospective-analysis' 

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
        newnc.variables[outtime].calendar = 'NOLEAP'
        #newnc.variables[outtime].cycle_length = np.float(365)
        newnc.variables[outtime][:] = dtot[nvar] + ndays[nmon]/2.0
        dtot[nvar]+=ndays[nmon]        

        newnc.createVariable(outvar[nvar], 'f', (outtime, 'lat', 'lon'), fill_value=spval[nvar])
        newnc.variables[outvar[nvar]].missing_value = spval[nvar]
        newnc.variables[outvar[nvar]].long_name = long_name[nvar]
        newnc.variables[outvar[nvar]].units = units[nvar]
        newnc.variables[outvar[nvar]].coordinates = 'lon lat'

        #get data
        var = fid.variables[invar[nvar]][:]
        if invar[nvar]== 'T2M':
           var+= -273.15 # Kelvin to Celcius
        newnc.variables[outvar[nvar]][:] = var
        newnc.close()
    fid.close()

