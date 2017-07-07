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
from lib_process_ERAinterim import ERAinterim_processing
import numpy as np
#import ConfigParser
import os 
import sys

# Define the list of interesting variables with their associated parameter on MARS server
ecmwf_param = {'u10' : '165.128', 'v10' : '166.128', 'd2' : '168.128' , 't2' : '167.128' , \
                'msl' : '151.128' , 'radsw' : '169.128' , 'radlw' : '175.128' , 'precip' : '228.128'}

# Choose years to download
fyear = 1980
lyear = 2010

server = ECMWFDataServer()

for year in np.arange(fyear,lyear+1):

	for key in ecmwf_param.keys():

		print 'working on variable', key, ' for year ', str(year)
		filegrb   =  key + '_ERAinterim_' + str(year) + '.grb'
                server.retrieve({
                    #'stream'    : "moda"
                    'stream'    : "oper", #oper = operational atmospheric model ( forecasting system used to generated the data when the same meteorological types are archived)
                    'levtype'   : "sfc",  #sfc = surface (denotes type of level)
                    'resol'     : "av",   #av = archived value for full spectral resolution
                    #'grid'      : "av",   #av = archived model grid
                    'grid'      : ".75/.75",
                    'param'     : ecmwf_param[key],
                    'dataset'   : "interim",
                    #'step'      : "0/3",
                    #'step'      : "3/6/9/12", ###
                    'time'      : "00/06/12/18", ###
                    'date'      : str(year) + "-01-01/to/" + str(year) + "-12-31",
                    'type'      : "an",
                    #'type'      : "fc",   #fc = forcast(type of fields to be retrieved)
                    'class'     : "ei",   #ei = ERA interim
                    'area'      : "55/-150/15/-105",     # Subset or clip to an area, here to California. Specify as North/West/South/East in Geographic lat/long degrees. Southern latitudes and western longitudes must be given as negative numbers.
                    #'format'    : "netcdf",
                    'target'    : filegrb 
                })
                #my_inputs = {'file':filenc,key:ecmwf_param[key],'year':year, 'ncumul':4, 'nx':512, 'ny':256}
                 #FIX NX, NY
                 #go = ERAinterim_processing(my_inputs)
                 #go()
