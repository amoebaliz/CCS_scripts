import numpy as np
import netCDF4 as nc

# OBJECTIVE: Getting ALL LENS endmemebers for SST surfPH and surfO2

# NCAR directory for Large Ensemble files
dir = '/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/'

# File name strings
his_fil_bas = 'b.e11.B20TRC5CNBDRD.f09_g16.'
fut_fil_bas = 'b.e11.BRCP85C5CNBDRD.f09_g16.'

# CCS DOMAIN
a = 245
b = 317
c = 225
d = 260

# Number of Years in timeseries
nyr = 2100-1950+1

# Number of NCAR LE endmembers
end_mem = 35

for VAR in ("SST", "O2", "PH"):
    fil_bas2 = '.pop.h.' + VAR + '.'
    # Initializing matrix for multi-endmember timeseries
    var_all_stor = np.empty((0,nyr*12,b-a,d-c))
    if VAR == "SST":
       endmembs = range(1,end_mem+1) 
    else: 
       endmembs = np.append((1,2),range(9,end_mem+1))
    for ne in endmembs:
        #######################
        ## HISTORICAL VALUES ##
        #######################
        if ne == 1: 
           his_str_yr = 1850
        else:
           his_str_yr = 1920

        # File name
        his_ncfil  = dir + VAR +'/' + his_fil_bas + str(ne).zfill(3) + fil_bas2 + str(his_str_yr) + '01-' + str(200512)+ '.nc'
        print his_ncfil
        fid_his = nc.Dataset(his_ncfil)

        #############
        # CESM GRID #
        #############
	if (ne == 1 and VAR == 'SST'):
           # SAVE LENS GRID INFO
           lat = fid_his.variables['TLAT'][a:b,c:d]
           lon = fid_his.variables['TLONG'][a:b,c:d]
           np.save('/glade/work/edrenkar/LENS_CCS_lat_lon.npy', np.append(np.expand_dims(lat,axis=0), np.expand_dims(lon,axis=0),axis=0))

        # Total number of yrs in historical file
        nyr_h = 2005-his_str_yr+1

        # From beginning of 1950 to end of length of timeseries
        Ifyr = 12*(nyr_h-(2005-1950+1))

        # Fetch Variables, PH lacks a depth dimension
        if VAR == 'PH':
           var_h = fid_his.variables[VAR][Ifyr:,a:b,c:d].squeeze()
        else:
           var_h = fid_his.variables[VAR][Ifyr:,0,a:b,c:d].squeeze()
        fid_his.close() 
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
            fut_ncfil = dir + VAR + '/' + fut_fil_bas + str(ne).zfill(3) + fil_bas2 + str(fut_str_yr[nf]) + '01-' + str(fut_end_yr[nf-nfils])+ '12.nc'
            fid_fut = nc.Dataset(fut_ncfil)
            # Number of years contained in future file 
            nyr_f = fut_end_yr[nf]-fut_str_yr[nf-nfils]+1
            # Monthly climatology/Annual Mean
            if VAR == 'PH':
               var_f = np.append(var_f, fid_fut.variables[VAR][:,a:b,c:d].squeeze(),axis=0)
            else:
               var_f = np.append(var_f, fid_fut.variables[VAR][:,0,a:b,c:d].squeeze(),axis=0)
        fid_fut.close()
        # STORE VALUES FOR ANALYSIS
        end_mem_vals = np.expand_dims(np.concatenate((var_h,var_f),axis=0),axis=0)
        var_all_stor = np.append(var_all_stor,end_mem_vals,axis=0)

    # SAVE FIELDS TO BINARY FILE
    npy_save_text = '/glade/work/edrenkar/LENS_1950-2100_' + VAR + '.npy'
    np.save(npy_save_text, var_all_stor)
