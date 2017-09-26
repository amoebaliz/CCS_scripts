from lib_process_ERAinterim import ERAinterim_processing
import calendar
import numpy as np
import sys
    
fyear = 1980
lyear = 2010

nday = [31,28,31,30,31,30,31,31,30,31,30,31]

dir_fc = '/Users/liz.drenkard/external_data/ERAinterim/Forecast/'
#dir_an = '/Users/liz.drenkard/external_data/ERAinterim/Analysis/'

for year in np.arange(fyear,lyear+1):
    if calendar.isleap(year):
       nday[1] = 29
    else:
       nday[1] = 28

    for mon in np.arange(1,12+1): 
        my_inputs = {'file_t2':     dir_fc + str(year) + '/t2_ERAinterim_'     + str(year) + '_' + str(mon).zfill(2) + '.nc.sub', \
                     'file_msl':    dir_fc + str(year) + '/msl_ERAinterim_'    + str(year) + '_' + str(mon).zfill(2) + '.nc.sub', \
                     'file_d2':     dir_fc + str(year) + '/d2_ERAinterim_'     + str(year) + '_' + str(mon).zfill(2) + '.nc.sub', \
                     'file_u10':    dir_fc + str(year) + '/u10_ERAinterim_'    + str(year) + '_' + str(mon).zfill(2) + '.nc.sub', \
                     'file_v10':    dir_fc + str(year) + '/v10_ERAinterim_'    + str(year) + '_' + str(mon).zfill(2) + '.nc.sub', \
                     'file_radlw':  dir_fc + str(year) + '/radlw_ERAinterim_'  + str(year) + '_' + str(mon).zfill(2) + '.nc.sub', \
                     'file_radsw':  dir_fc + str(year) + '/radsw_ERAinterim_'  + str(year) + '_' + str(mon).zfill(2) + '.nc.sub', \
                     'file_precip': dir_fc + str(year) + '/precip_ERAinterim_' + str(year) + '_' + str(mon).zfill(2) + '.nc.sub', \
                     'year':year, 'month':str(mon).zfill(2), 'ndays':nday[mon-1], 'nx':80, 'ny':80, 'output_dir':dir_fc}

        go = ERAinterim_processing(my_inputs)
        go()
