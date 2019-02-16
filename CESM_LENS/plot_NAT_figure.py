import numpy as np
import netCDF4 as nc

# OBJECTIVE: Getting ALL LENS endmemebers for SST surfPH and surfO2

# NCAR directory for Large Ensemble files
dir = '/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/'

# File name strings
his_fil_bas = 'b.e11.B20TRC5CNBDRD.f09_g16.'
fut_fil_bas = 'b.e11.BRCP85C5CNBDRD.f09_g16.'

fil_bas2 = '.pop.h.SST.'

# CCS DOMAIN
a = 245
b = 317
c = 225
d = 260

# Number of Years in timeseries
nyr = 2100-1950+1

# Number of NCAR LE endmembers
end_mem = 35

for VAR in ("SST", "PH", "O2"):
    # Initializing matrix for multi-endmember timeseries
    var_stor_all=np.empty((0,ny*12,b-a,d-c))

    for ne in range(1,end_mem+1):
        #######################
        ## HISTORICAL VALUES ##
        #######################
        if ne == 1: 
           his_str_yr = 1850
        else:
           his_str_yr = 1920

        # File name
        his_ncfil  = dir + his_fil_bas + str(ne).zfill(3) + fil_bas2 + str(his_str_yr) + '01-' + str(200512)+ '.nc'
        fid_his = nc.Dataset(his_ncfil)

        #############
        # CESM GRID #
        #############
	if (ne == 1 and VAR == 'SST'):
           # SAVE LENS GRID INFO
           lat = fid_his.variables['TLAT'][a:b,c:d]
           lon = fid_his.variables['TLONG'][a:b,c:d]
           np.save('LENS_CCS_lat_lon.npy', np.ma.append(np.ma.expand_dims(lat,axis=0), np.ma.expand_dims(lon,axis=0),axis=0))

        # Total number of yrs in historical file
        nyr_h = 2005-his_str_yr+1

        # From beginning of 1950 to end of length of timeseries
        Ifyr = 12*(nyr_h-(2005-1950+1))

        # Fetch Variables, PH lacks a depth dimension
        if VAR == 'PH':
           var_h = fid_his.variables[VAR][Ifyr:,a:b,c:d].squeeze()
        else:
           var_h = fid_his.variables[VAR][Ifyr:,0,a:b,c:d].squeeze()
        
        ################### 
        ## FUTURE VALUES ##
        ###################
        fut_str_yr = [2006,2081]
        fut_end_yr = [2080,2100]

        if ne < 34:
           # Two files: 2006-2080, 2081-2100
           tot_yr_f = [np.diff(fut_end_yr),np.diff(fut_end_yr)]
        else:
           # One file: 2006-2100
           tot_yr_f = [(2100-2006) + 1]

        # Opening and appending variable number of files
        var_f = np.empty((0,b-a,d-c))
        nfils = len(tot_yr_f)
        for nf in range(nfils):
            fut_ncfil = dir + fut_fil_bas + str(ne).zfill(3) + fil_bas2 + str(fut_str_yr[nf]) + '01-' + str(fut_end_yr[nf-nfils])+ '12.nc'
            fid_fut = nc.Dataset(fut_ncfil)
            # Number of years contained in future file 
            nyr_f = fut_end_yr[nf]-fut_str_yr[nf-nfils]+1
            # Monthly climatology/Annual Mean
            if VAR == 'PH':
               var_f = np.ma.append(var_f, fid_fut.variables[VAR][:,a:b,c:d].squeeze())
            else:
               var_f = np.ma.append(var_f, fid_fut.variables[VAR][:,0,a:b,c:d].squeeze())

        # STORE VALUES FOR ANALYSIS 
        var_all_stor = np.ma.append(var_all_stor,np.ma.expand_dims(np.ma.append(var_h,var_f),axis=0),axis=0)

    # SAVE FIELDS TO BINARY FILE
    npy_save_text = 'LENS_' + VAR + '.npy'
    np.save(npy_save_text, var_all_stor)
