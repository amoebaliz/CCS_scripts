import numpy as np
import netCDF4 as nc
import csv
import struct


grdfil = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid = nc.Dataset(grdfil)

lon_rho  = fid.variables['lon_rho'][:]
lat_rho  = fid.variables['lat_rho'][:]
mask      = fid.variables['mask_rho'][:]
h         = fid.variables['h'][:] 

print len(lat_rho)

sfname = 'California_Current_System_Domain_Bathymetry'
with open(sfname, 'w') as f:
        writer = csv.writer(f, delimiter=' ', lineterminator='\n')
        for jj in range(len(lon_rho)):
            row = h[jj,:]
            writer.writerow(row)

sfname = 'California_Current_System_Domain_Latitude'
with open(sfname, 'w') as f:
        writer = csv.writer(f, delimiter=' ', lineterminator='\n')
        for jj in range(len(lon_rho)):
            row = lat_rho[jj,:]
            writer.writerow(row)

sfname = 'California_Current_System_Domain_Longitude'
with open(sfname, 'w') as f:
        writer = csv.writer(f, delimiter=' ', lineterminator='\n')
        for jj in range(len(lon_rho)):
            row = lon_rho[jj,:]
            writer.writerow(row)

sfname = 'California_Current_System_Domain_Mask'
with open(sfname, 'w') as f:
        writer = csv.writer(f, delimiter=' ', lineterminator='\n')
        for jj in range(len(lon_rho)):
            row = mask[jj,:]
            writer.writerow(row)



