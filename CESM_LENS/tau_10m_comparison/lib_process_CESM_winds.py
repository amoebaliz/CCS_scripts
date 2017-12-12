import netCDF4 as nc
import numpy as np
import calendar
import datetime as dt
import wind_toolbox
import os

def write_ncfile(lon_array,lat_array,time,var,dict_wrt):
                fid = nc.Dataset(dict_wrt['fileout'], 'w', format='NETCDF3_CLASSIC')
                fid.description = 'CESM post-processing based on code by raphael.dussin@gmail.com'
                # dimensions
                fid.createDimension('lat', lat_array.shape[0])
                fid.createDimension('lon', lon_array.shape[0])
                fid.createDimension(dict_wrt['time_dim'], None)
                # variables
                latitudes  = fid.createVariable('lat', 'f8', ('lat',))
                longitudes = fid.createVariable('lon', 'f8', ('lon',))
                times      = fid.createVariable(dict_wrt['time_var'], 'f8', (dict_wrt['time_dim'],))
                variable   = fid.createVariable(dict_wrt['varname'], 'f4', (dict_wrt['time_dim'],'lat','lon',),fill_value=1e15)

                # attributes
                longitudes.units = "degrees_east"
                longitudes.valid_min = lon_array.min()
                longitudes.valid_max = lon_array.max()
                longitudes.long_name = "longitude"

                latitudes.units = "degrees_north"
                latitudes.valid_min = lat_array.min()
                latitudes.valid_max = lat_array.max()
                latitudes.long_name = "latitude"

                times.units = "days" 
                #times.valid_min = time.min()
                #times.valid_max = time.max()
                times.calendar = "NO LEAP"

                variable.long_name = dict_wrt['long name']
                variable.units = dict_wrt['units']
                variable.coordinates = "lon lat"
                variable.time = dict_wrt['time_var']
                variable.missing_value = 1e15
                #variable.valid_range = var.min() , var.max()

                #data
                latitudes[:]    = lat_array
                longitudes[:]   = lon_array
                times[:]        = time
                variable[:,:,:] = var

                # close
                fid.close()
                return None

#dir = '/glade/p/work/edrenkar/external_data/LENS/scripts/his/'
dir = '/glade/p/work/edrenkar/external_data/LENS/scripts/016/'

# HISTORICAL
tx_fil  = dir + 'TAUX_81-05_clim.nc'
ty_fil  = dir + 'TAUY_81-05_clim.nc'
u10_fil = dir + 'U10_81-05_clim.nc'

uout_fil = dir + 'Uwind_81-05_clim.nc'
vout_fil = dir + 'Vwind_81-05_clim.nc'

print uout_fil
print vout_fil

# FUTURE
tx_fil  = dir + 'TAUX_RCP85_016_76-00_clim.nc' 
ty_fil  = dir + 'TAUY_RCP85_016_76-00_clim.nc'
u10_fil = dir + 'U10_RCP85_016_76-00_clim.nc' 

uout_fil = dir + 'Uwind_016_76-00_clim.nc' 
vout_fil = dir + 'Vwind_016_76-00_clim.nc'

print uout_fil
print vout_fil

fidtx = nc.Dataset(tx_fil)
fidty = nc.Dataset(ty_fil)
fidu10 = nc.Dataset(u10_fil)

taux = fidtx.variables['TAUX'][:].squeeze()
tauy = fidty.variables['TAUY'][:].squeeze()
u10  = fidu10.variables['U10'][:].squeeze()

lat = fidu10.variables['lat'][:]
lon = fidu10.variables['lon'][:]
time = fidu10.variables['time'][:]

u_out = np.zeros(u10.shape)
v_out = np.zeros(u10.shape)

for nt in range(12):
    u_tmp, v_tmp = wind_toolbox.u_v_from_taux_tauy_u10(taux[nt,:].squeeze(),tauy[nt,:].squeeze(),u10[nt,:].squeeze())

    u_out[nt,:] = -1*u_tmp 
    v_out[nt,:] = -1*v_tmp

my_dictu = {'varname':'Uwind','time_dim':'time','time_var':'time','long name':'Zonal Velocity at 10m',\
                        'units':'m s-1','fileout':uout_fil}

my_dictv = {'varname':'Vwind','time_dim':'time','time_var':'time','long name':'Meridional Velocity at 10m',\
                        'units':'m s-1','fileout':vout_fil}        

write_ncfile(lon,lat[::-1],time,u_out[:,::-1,:],my_dictu)
write_ncfile(lon,lat[::-1],time,v_out[:,::-1,:],my_dictv)
