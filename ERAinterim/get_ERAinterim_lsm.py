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
#from lib_process_ERAinterim import ERAinterim_processing
import calendar
import numpy as np
import os 
import sys

server = ECMWFDataServer()
grbfil = 'ERAinterim_lsm.grb' 
ncfil  = 'ERAinterim_lsm.nc'
          
# MARS retreival request
server.retrieve({
                    'stream'    : "oper", #oper = operational atmospheric model ( forecasting system used to generated the data when the same meteorological types are archived)
                    'levtype'   : "sfc",  #sfc = surface (denotes type of level)
                    'resol'     : "av",   #av = archived value for full spectral resolution
                    'param'     : "172.128",
                    'dataset'   : "interim",
                    'step'      : "0",
                    #'time'      : "00/06/12/18", ###
                    'time'      : "1200", ###
                    'date'      : "19890101",
                    'type'      : "an",   #fc = forcast(type of fields to be retrieved)
                    'class'     : "ei",   #ei = ERA interim
                    'target'    : grbfil
                })

                
# CONVERT grib to nc file
grb2nc = 'cdo -R -t ecmwf -f nc -r copy ' + grbfil + ' ' + ncfil               
os.system(grb2nc)
cln_up = 'rm ' + grbfil
os.system(cln_up)
                
                
# SPATIALLY subset nc file to CCS
CCS_sub = 'ncks -O -d lat,40,119 -d lon,290,369 ' + ncfil + ' /Users/elizabethdrenkard/external_data/' + ncfil + '.sub'
os.system(CCS_sub)
cln_up = 'rm ' + ncfil 
os.system(cln_up)

                # PROCESS ncfil to be ROMS compatible
               
                 #my_inputs = {'file':filenc,key:ecmwf_param[key],'year':year, 'ncumul':4, 'nx':80, 'ny':80}
                 #go = ERAinterim_processing(my_inputs)
                 #go()
