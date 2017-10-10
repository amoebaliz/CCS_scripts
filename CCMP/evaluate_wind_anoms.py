import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pyroms

def make_mask():
    lat = fid.variables['latitude'][:]
    lon = fid.variables['longitude'][:]
    mask = np.zeros((len(lat),len(lon)))

    # CCS domain corners
    # SW: ccslat[ 0, 0]; ccslon[ 0, 0]
    # SE: ccslat[ 0,-1]; ccslon[ 0,-1]
    # NW: ccslat[-1, 0]; ccslon[-1, 0]
    # NE: ccslat[-1,-1]; ccslon[-1,-1]

    # polynomial equations
    p_gt = np.poly1d(np.polyfit(ccslat[:,0],ccslon[:,0],2))
    p_lt = np.poly1d(np.polyfit(ccslat[0,:],ccslon[0,:],2))

   # Iterate over all latitudes included in domain
    for jj in range(len(lat)):
        # if pass SE and NW corners, update indices for slopes & intercepts  
        if lat[jj] > ccslat[0,-1]:
           p_lt = np.poly1d(np.polyfit(ccslat[:,-1],ccslon[:,-1],2))
           if lat[jj] > ccslat[-1,0]: 
              p_gt = np.poly1d(np.polyfit(ccslat[-1,:],ccslon[-1,:],2)) 

        # if lon is greater than gt_val and less than lt_val, mask = 1
        ilon = np.where((p_gt(lat[jj]) < lon) & (lon < p_lt(lat[jj])))[0]
        mask[jj,ilon]=1

    # fill in border points to account for non-linearity of ROMS grid
    for nt in range(1):
        mask_sub = np.zeros(mask.shape)
        for jj in range(1,len(lat)-1):
            for ji in range(1,len(lon)-1):
                # Most needed on souther and western boarders
                # :. Add more only if there are masked points above and/or to the right
                mask_sum = np.sum([mask[jj,ji+1],mask[jj+1,ji]])
                if mask[jj,ji]==0 and mask_sum>0:
                   mask_sub[jj,ji]=1  
        mask[mask_sub>0]=1

    return lat, lon, mask
########################################

grd = pyroms.grid.get_ROMS_hgrid('CCS')
ccslat = grd.lat_rho
ccslon = grd.lon_rho

lat_max=np.max(ccslat)
lat_min=np.min(ccslat)

n = 0
# load individual CCSM files
dir = '/glade/p/work/edrenkar/external_data/CCMP/'
ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
# LOOP OVER ALL YEARS
for year in range(1990,2010+1):
    print year
    # LEAP YEAR EVALUATION
    ineq = year%4
    if ineq == 0:
       ndays[1] = 29
    else:
       ndays[1] = 28

    # LOOP OVER ALL MONTHS
    for mon in range(12):
        # LOOP OVER ALL DAYS
        for day in range(ndays[mon]):
            ncfil= dir + str(year)+ '/CCMP_Wind_Analysis_'+ str(year) + str(mon+1).zfill(2)+str(day+1).zfill(2) + '_daily_anom.nc'
            print ncfil
            # first time, load/store the CCSM lat lon
            fid = nc.Dataset(ncfil)
            uanom = fid.variables['u_anom'][:].squeeze()
            vanom = fid.variables['v_anom'][:].squeeze()
            if n < 1:
               # Build a mask of 1's falling withing the CCD grid domain
               lat, lon, ccmp_mask = make_mask()
               ustore = np.array((year, mon+1, day+1,\
                        np.mean(abs(uanom[ccmp_mask>0])),\
                        np.std(abs(uanom[ccmp_mask>0])), \
                        np.min(abs(uanom[ccmp_mask>0])), \
                        np.max(abs(uanom[ccmp_mask>0]))))

               vstore = np.array((year, mon+1, day+1,\
                        np.mean(abs(vanom[ccmp_mask>0])),\
                        np.std(abs(vanom[ccmp_mask>0])), \
                        np.min(abs(vanom[ccmp_mask>0])), \
                        np.max(abs(vanom[ccmp_mask>0]))))
               n+=1

            bu = np.array((year, mon+1, day+1,\
                           np.mean(abs(uanom[ccmp_mask>0])),\
                           np.std(abs(uanom[ccmp_mask>0])), \
                           np.min(abs(uanom[ccmp_mask>0])), \
                           np.max(abs(uanom[ccmp_mask>0]))))

            bv = np.array((year, mon+1, day+1,\
                           np.mean(abs(vanom[ccmp_mask>0])),\
                           np.std(abs(vanom[ccmp_mask>0])), \
                           np.min(abs(vanom[ccmp_mask>0])), \
                           np.max(abs(vanom[ccmp_mask>0]))))
 
            ustore = np.concatenate((ustore,bu),axis=0)
            vstore = np.concatenate((vstore,bv),axis=0) 

np.savetxt('uwind_anomaly_stats', ustore, delimiter=',') 
np.savetxt('vwind_anomaly_stats', vstore, delimiter=',') 
