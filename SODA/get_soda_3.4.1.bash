#!/bin/bash
#
while read f; do
    wget $f
    ncfil=${f#*ocean/}
    echo $ncfil
    ncks -O -d xt_ocean,501,700 -d yt_ocean,541,740 -d xu_ocean,501,700 -d yu_ocean,541,740 $ncfil ../${ncfil}.sub
    rm $ncfil
done <soda3.4.1_5dy_ocean_or_original.txt 
