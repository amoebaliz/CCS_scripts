import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import ESMF

# ERAint files
# yellowstone
#fidUclim = nc.Dataset('/glade2/scratch2/edrenkar/CCS-inputs/drowned_ERAi_u10_1981-2010_monthly_clim.nc')
#fidVclim = nc.Dataset('/glade2/scratch2/edrenkar/CCS-inputs/drowned_ERAi_v10_1981-2010_monthly_clim.nc')

# variable file names
avars = ['msl','q2','radlw','radsw','t2']
#avars = ['precip']
Var_nm = ['Pair','Qair','lwrad_down','swrad','Tair']
#Var_nm = ['rain']

lens_dir = '/Users/elizabethdrenkard/Desktop/016/'
ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
time_vals = np.zeros(12)
dtot = 0 

for n in range(12):
    time_vals[n] = ndays[n]/2.0 + dtot
    dtot += ndays[n]

for nv in range(len(avars)):
    print Var_nm[nv]
    Encfil = lens_dir + 'drowned_ERAi_' + avars[nv] + '_1981-2010_monthly_clim.nc'
    fid_clim = nc.Dataset(Encfil)
    era_var = fid_clim.variables[Var_nm[nv]][:].squeeze()

    Lncfil = lens_dir + 'drowned_LENS_016_' + avars[nv] + '_delta.nc'
    fid_delta = nc.Dataset(Lncfil)

    if nv == 2:
       Var_nm[nv] = 'lwrad'

    lens_var =  fid_delta.variables[Var_nm[nv]][:].squeeze()

    if nv == 2:
       Var_nm[nv] = 'lwrad_down'

    dlon = fid_delta.variables['lon'][:]
    dlon[dlon>180]=dlon[dlon>180]-360
    dlat = fid_delta.variables['lat'][:]

    Xi,Yi = np.meshgrid(dlon,dlat)
    sourcegrid = ESMF.Grid(np.array(Xi.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)
    #sourcegrid = ESMF.Grid(np.array(dlon.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)


    lon = fid_clim.variables['lon'][:]
    lon[lon>180]=lon[lon>180]-360
    lat = fid_clim.variables['lat'][:]

    nx = len(lon)
    ny = len(lat)

    Xn, Yn = np.meshgrid(lon,lat)
    destgrid = ESMF.Grid(np.array(Xn.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

    ## POINTERS
    source_lon = sourcegrid.get_coords(0)
    source_lat = sourcegrid.get_coords(1)
    dest_lon = destgrid.get_coords(0)
    dest_lat = destgrid.get_coords(1)

    ## FILLS
    source_lon[...] = Xi
    source_lat[...] = Yi
    #source_lon[...] = dlon # precip
    #source_lat[...] = dlat # precip
    dest_lon[...] = Xn
    dest_lat[...] = Yn

    sourcefield = ESMF.Field(sourcegrid, name = 'CESM_Delta')
    destfield = ESMF.Field(destgrid, name = 'ERAi_Clim')

    # Allocate output variables
    var_out = np.zeros(era_var.shape)

    ## ITERATE OVER ALL MONTHS 
    for nmon in range(12):
        # REGULAR GRID BILINIEAR INTERPOLATION
        sourcefield.data[...] = lens_var[nmon,:].squeeze()
        regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.BILINEAR,  
                     unmapped_action = ESMF.UnmappedAction.IGNORE)
        destfield = regrid(sourcefield, destfield) 
   
        var_out[nmon,:] = era_var[nmon,:].squeeze() + destfield.data
 
    # Save new wind files
    ncfile = 'ERAi_CESM_016_delta_' + avars[nv] + '.nc'
    fid2 = nc.Dataset(ncfile,'w')

    fid2.createDimension('time', None)
    fid2.createDimension('lat', ny)
    fid2.createDimension('lon', nx)

    fid2.createVariable('time', 'f8', ('time'))
    fid2.variables['time'].units = fid_clim.variables['time'].units
    fid2.variables['time'].cycle_length = fid_clim.variables['time'].cycle_length
    fid2.variables['time'][:] = time_vals

    fid2.createVariable('lat','f8',('lat'))
    fid2.variables['lat'].long_name = fid_clim.variables['lat'].long_name
    fid2.variables['lat'].units = fid_clim.variables['lat'].units
    fid2.variables['lat'][:]=lat

    fid2.createVariable('lon','f8',('lon'))
    fid2.variables['lon'].long_name = fid_clim.variables['lon'].long_name
    fid2.variables['lon'].units = fid_clim.variables['lon'].units
    fid2.variables['lon'][:]=lon
    
    fid2.createVariable(Var_nm[nv],'f8',('time','lat','lon'),fill_value = np.float(1.0e15))
    fid2.variables[Var_nm[nv]].long_name = fid_clim.variables[Var_nm[nv]].long_name
    fid2.variables[Var_nm[nv]].units = fid_clim.variables[Var_nm[nv]].units
    fid2.variables[Var_nm[nv]].coordinates = fid_clim.variables[Var_nm[nv]].coordinates
    fid2.variables[Var_nm[nv]].time = fid_clim.variables[Var_nm[nv]].time
    u_txt = "ERAinterim climatology (1981-2010) + CESM Delta"
    fid2.variables[Var_nm[nv]].details = u_txt
    fid2.variables[Var_nm[nv]][:]=var_out

    fid2.close()
    fid_clim.close()
