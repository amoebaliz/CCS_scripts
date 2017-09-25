#!/usr/bin/env python
#
# (C) Copyright 2012-2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.
#
# To run this example, you need an API key 
# available from https://api.ecmwf.int/v1/key/
#
# R.Dussin 2014 - adapt script to download ERAinterim
# L.Drenkard 2017 - additional adaptations

from ecmwfapi import ECMWFDataServer
import calendar
import numpy as np
import os 

# Define the list of interesting variables with their associated parameter on MARS server
#ecmwf_param = {'u10' : '165.128', 'v10' : '166.128', 'd2' : '168.128' , 't2' : '167.128' , \
#                'msl' : '151.128' , 'radsw' : '169.128' , 'radlw' : '175.128' , 'precip' : '228.128'}

# NOTE: need to use net short wave radiation (not downward) as above - different code
ecmwf_param = {'radsw' : '176.128'}

# Choose years to download
fyear = 1981
lyear = 2010
nday = [31,28,31,30,31,30,31,31,30,31,30,31]

server = ECMWFDataServer()

for year in np.arange(fyear,lyear+1):
    if calendar.isleap(year):
       nday[1] = 29
    else: 
       nday[1] = 28 
    for mon in np.arange(1,12+1):
	for key in ecmwf_param.keys():

		print 'working on variable', key, ' for year ', str(year)
		grbfil   =  key + '_ERAinterim_' + str(year) + '_' + str(mon).zfill(2) + '.grb'
                ncfil =  key + '_ERAinterim_' + str(year) + '_' + str(mon).zfill(2) + '.nc'

                # MARS retreival request
                server.retrieve({
                    'stream'    : "oper", #oper = operational atmospheric model ( forecasting system used to generated the data when the same meteorological types are archived)
                    'levtype'   : "sfc",  #sfc = surface (denotes type of level)
                    'resol'     : "av",   #av = archived value for full spectral resolution
                    'param'     : ecmwf_param[key],
                    'dataset'   : "interim",
                    'step'      : "3/6/9/12", ###
                    #'time'      : "00/06/12/18", ###
                    'time'      : "00/12", ###
                    'date'      : str(year) + "-" + str(mon).zfill(2) + "-01/to/" + str(year) + "-" + str(mon).zfill(2) + "-" + str(nday[mon-1]).zfill(2),
                    #'type'      : "an",
                    'type'      : "fc",   #fc = forcast(type of fields to be retrieved)
                    'class'     : "ei",   #ei = ERA interim
                    'target'    : grbfil 
                })

                # CONVERT grib to nc file
                grb2nc = 'cdo -R -t ecmwf -f nc -r copy ' + grbfil + ' ' + ncfil               
                os.system(grb2nc)
                cln_up = 'rm ' + grbfil
                os.system(cln_up)
                
                # SPATIALLY subset nc file to CCS
                CCS_sub = 'ncks -O -d lat,40,119 -d lon,290,369 ' + ncfil + ' ' + ncfil + '.sub'
                os.system(CCS_sub)
                cln_up = 'rm ' + ncfil 
                os.system(cln_up)

                # PROCESS ncfil to be ROMS compatible
               
                 #my_inputs = {'file':filenc,key:ecmwf_param[key],'year':year, 'ncumul':4, 'nx':80, 'ny':80}
                 #go = ERAinterim_processing(my_inputs)
                 #go()
