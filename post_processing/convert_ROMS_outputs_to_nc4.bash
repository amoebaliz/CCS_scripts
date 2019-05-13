#!/bin/bash

# NOTHING TO TOUCH HERE !!!

# check input file
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 rc-file" >&2
  echo "Example: $0 example-CCS1-cobalt7" >&2
  echo "Provide all necessary informations for your run in rcfile and run me again"
  exit 1
fi
if ! [ -f "$1" ]; then
  echo "rc-file $1 not found" >&2
  exit 1
fi

# Load informations about run
rcfile=$1
. ./$rcfile

if [[ -z $DIRROOT ]] ;  then echo no valid DIRROOT provided, aborting ; exit 1 ; fi
if [[ -z $RUN ]] ;      then echo no valid RUN provided, aborting     ; exit 1 ; fi
if [[ -z $FYEAR ]] ;    then echo no valid FYEAR provided, aborting   ; exit 1 ; fi
if [[ -z $LYEAR ]] ;    then echo no valid LYEAR provided, aborting   ; exit 1 ; fi

if [ ! -d $DIRROOT/$RUN ] ; then echo $DIRROOT/$RUN does not exist ; exit 1 ; fi

cd $DIRROOT/$RUN/

# loop on year directories
for yr in $( seq $FYEAR $LYEAR ) ; do
    year=$(printf %04d $yr)
    if [ ! -d $DIRROOT/$RUN/$year ] ; then echo $DIRROOT/$RUN/$year does not exist ; exit 1 ; fi
    cd $DIRROOT/$RUN/$year

    # loop on files
    for file in $( ls *.nc ) ; do

        echo working on file $file from year $year
        mv $file ${file}3
        ncks -4 -L 1 ${file}3 -o $file
        if [ -f $file ] ; then
           ncatted -O -a format,global,m,c,'netCDF-4 compression level 1' $file
           rm ${file}3
        else
           echo 'netcdf conversion failed for file' $file
           mv ${file}3 $file
        fi

    done

done
