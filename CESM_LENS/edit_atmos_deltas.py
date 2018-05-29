import numpy as np
import netCDF4 as nc

# MAKING VARIABLES AND FILENAMES SIMILAR TO ERAinterim


def build_atmos_fil(run,nvar):
    runstr = str(run).zfill(3)

    # OLD CLIM FILE 
    infil = dir + runstr + '_' + vars[nvar] + '_clim_delta.nc'      
    fidin = nc.Dataset(infil)
    TIME = fidin.variables['time']
    #LAT = fidin.variables['TLAT']
    #LON = fidin.variables['TLONG']
    LAT = fidin.variables['lat']
    LON = fidin.variables['lon']
    OLD_VAR = fidin.variables[vars[nvar]]

    # NEW FILE BEING CREATED
    ncfil = '/glade/p/work/edrenkar/Inputs/construct/Deltas/LENS_deltas_atm/' \
            + runstr + '/' + 'LENS_' + runstr +'_' + tit_vars[nvar] + '_delta.nc'
    fid = nc.Dataset(ncfil,'w')

    # Dimensions 
    fid.createDimension('time', None)
    fid.createDimension('lat',50)
    fid.createDimension('lon',50)
    #fid.createDimension('lat',80)
    #fid.createDimension('lon',80)

    # TIME
    fid.createVariable('time', 'f8', ('time'))
    fid.variables['time'][:] = TIME[:]
      
    # LAT
    fid.createVariable('lat', 'f8', ('lat'))
    #fid.createVariable('lat', 'f8', ('lat','lon'))
    fid.variables['lat'].long_name = LAT.long_name 
    fid.variables['lat'].units = LAT.units   
    fid.variables['lat'][:] = LAT[:]

    # LON
    fid.createVariable('lon', 'f8', ('lon')) 
    #fid.createVariable('lon', 'f8', ('lat','lon'))  
    fid.variables['lon'].long_name = LON.long_name
    fid.variables['lon'].units = LON.units
    fid.variables['lon'][:] = LON[:]

    # VARIABLE
    fid.createVariable(new_vars[nvar], 'f8', ('time','lat','lon'))
    fid.variables[new_vars[nvar]].units = OLD_VAR.units
    fid.variables[new_vars[nvar]].long_name = OLD_VAR.long_name
    fid.variables[new_vars[nvar]][:] = OLD_VAR[:].squeeze()

    fid.close() 
    print ncfil
# -----------
vars= ['Q','PSL','TREFHT','U','V','FSDS','FLDS'] #,'PREC_F'}
#vars = ['PREC_F']
tit_vars = ['q2','msl','t2','u10','v10','radsw','radlw'] #,'precip'}
#tit_vars = ['precip']
new_vars = ['Qair','Pair','Tair','Uwind','Vwind','swrad', 'lwrad'] #,'rain')
#new_vars = ['rain']

dir = '/glade/p/work/edrenkar/external_data/LENS/difs/'

for run in [17]:
    runstr = str(run).zfill(3)
    for nvar in range(len(vars)):
        # CREATE NEW MONTHLY CLIM FILE WITH ALL OCEAN FIELDS
        build_atmos_fil(run,nvar)


