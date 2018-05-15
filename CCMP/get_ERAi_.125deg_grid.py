#!/usr/bin/env python
import os  
import sys 
from ecmwfapi import ECMWFDataServer

grbfil =  'ERAinterim_.125_lsm.grb'
ncfil  =  'ERAinterim_.125_lsm.nc'

server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "1989-01-01",
    "expver": "1",
    "grid": "0.125/0.125",
    "levtype": "sfc",
    "param": "172.128",
    "step": "0",
    "stream": "oper",
    "time": "12:00:00",
    "type": "an",
    "target": grbfil,
})

# CONVERT grib to nc file
grb2nc = 'cdo -R -t ecmwf -f nc -r copy ' + grbfil + ' ' + ncfil    
os.system(grb2nc)
# REMOVE grib file
cln_up = 'rm ' + grbfil
os.system(cln_up)
