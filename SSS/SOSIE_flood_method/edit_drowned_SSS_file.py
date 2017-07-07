import netCDF4 as nc
from datetime import datetime 

ncfil1 = 'sss_monthly_climatology.nc'
ncfil2 = 'drowned_sss_monthly_climatology.nc'
    
fid1 = nc.Dataset(ncfil1)
fid2 = nc.Dataset(ncfil2,'a')

fid2.variables['time'].units = fid1.variables['sss_time'].units
#fid2.variables['time'].cycle_length = fid1.variables['sss_time'].cycle_length

fid2.variables['lon'].long_name = fid1.variables['lon'].long_name
fid2.variables['lat'].long_name = fid1.variables['lat'].long_name

fid2.variables['SSS'].missing_value = fid1.variables['SSS'].missing_value

fid1.close()

#create ROMS forcing file
fid2.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
fid2.conventions = 'CF-1.6'
fid2.title = 'A Global Ocean Hydrography with a High Quality Arctic Ocean Version 3.0'
fid2.reference = 'Steele, Morley and Ermold (2001), PHC: A global ocean hydrography with a high-quality Arctic Ocean, J. Clim'
fid2.source = 'Polar Science Center/Applied Physics Lab/University of Washington'
fid2.institution = 'University of Washington'

fid2.close()

