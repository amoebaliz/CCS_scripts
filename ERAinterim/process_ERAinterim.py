from lib_process_ERAinterim import ERAinterim_processing
import calendar
import numpy as np
import sys
import os   

fyear = 1982
lyear = 1982

nday = [31,28,31,30,31,30,31,31,30,31,30,31]
#MACBOOK PRO
dir_fc = '/Users/elizabethdrenkard/ANALYSES/CCS/ERAinterim/'
# NOAA SWFSC
dir_fc = '/Users/liz.drenkard/external_data/ERAinterim/Forecast/'
for year in np.arange(fyear,lyear+1):
    if calendar.isleap(year):
       nday[1] = 29
    else:
       nday[1] = 28

    for mon in np.arange(1,1+1): 
        my_inputs = {'file_t2':     't2_ERAinterim_'     + str(year) + '_' + str(mon).zfill(2) + '.nc', \
                     'file_msl':    'msl_ERAinterim_'    + str(year) + '_' + str(mon).zfill(2) + '.nc', \
                     'file_d2':     'd2_ERAinterim_'     + str(year) + '_' + str(mon).zfill(2) + '.nc', \
                     'file_u10':    'u10_ERAinterim_'    + str(year) + '_' + str(mon).zfill(2) + '.nc', \
                     'file_v10':    'v10_ERAinterim_'    + str(year) + '_' + str(mon).zfill(2) + '.nc', \
                     'file_radlw':  'radlw_ERAinterim_'  + str(year) + '_' + str(mon).zfill(2) + '.nc', \
                     'file_radsw':  'radsw_ERAinterim_'  + str(year) + '_' + str(mon).zfill(2) + '.nc', \
                     'file_precip': 'precip_ERAinterim_' + str(year) + '_' + str(mon).zfill(2) + '.nc', \
                     'year':year, 'month':str(mon).zfill(2), 'ndays':nday[mon-1], 'ncumul':4,'nx':512, 'ny':256, 'output_dir':dir_fc}

        #my_inputs = {'file_radsw':  dir_fc + '/radsw_ERAinterim_' + str(year) + '_' + str(mon).zfill(2) + '.nc.sub',\
        #             'year':year, 'month':str(mon).zfill(2), 'ndays':nday[mon-1], 'ncumul':4,'nx':80, 'ny':80, 'output_dir':dir_fc}


        go = ERAinterim_processing(my_inputs)
        go()

        cln_up = 'rm ' + ncfil 
        os.system(cln_up)
