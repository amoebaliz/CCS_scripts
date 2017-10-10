#!/bin/bash

for year in $( seq 1980 2010 ) ; do

    mkdir ../$year
    mv soda3.3.1_mon_${year}*nc ./$year/.
    #mv ./$year/* .
done

