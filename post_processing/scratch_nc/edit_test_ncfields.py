import numpy as np
import netCDF4 as nc

nc_vars = ['v','u']

for nt in range(4):
    ncfil = 'test_' + str(nt+1).zfill(3) + '.nc'
    fid = nc.Dataset(ncfil,'a')

    fid.variables['ocean_time'][:]=nt

    for var in nc_vars:
        fid.variables[var][:] = nt + 2
        #fid.variables[var][:] = np.random.randn(var_ob.shape[0],var_ob.shape[1])

fid.close()



