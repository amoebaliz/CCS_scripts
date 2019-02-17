import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import datetime as dt
# OBJECTIVE: Plot Surf Temp, O2, PH time-series
dir = '/glade/work/edrenkar/'
grid_fil = dir + 'LENS_CCS_lat_lon.npy'

latlon = np.load(grid_fil)
lat = latlon[0,:].squeeze()
lon = latlon[1,:].squeeze()

# Number of Years in timeseries
nyrs = 2100-1950+1

# ANNUAL or MONTHLY MEANS
# ANNUAL = 0, MONTHLY = 1
an_val = 0

# Assign Month-Center Dates
ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
if an_val:
   cesm_dates = np.zeros(nyrs*12)
else:
   cesm_dates = np.zeros(nyrs)   
n=0

for nyr in range(1950,2100+1):

    if an_val:
       if nyr%4 == 0:
          ndays[1] = 29
       else:
          ndays[1] = 28
       for nmon in range(12):
          date_val = dt.datetime(nyr,nmon+1,1)
          cesm_dates[n] = pltd.date2num(date_val + dt.timedelta(days=ndays[nmon]/2.-1))
          n+=1
    else:
       cesm_dates[n] = pltd.date2num(dt.datetime(nyr,1,1) + dt.timedelta(days=365/2.-1))
       n+=1 

# Figure setup:SST, PH, O2
fig, ax = plt.subplots(3, figsize=(5,8),sharex=True)
ytic = [range(17,23+1,3), np.arange(7.7,8.1+.1,.2), range(229,255+1,11)]        
nv=0
for VAR in ("SST", "PH", "O2"):
    npyfil = dir + 'LENS_1950-2100_' + VAR + '.npy'
    var = np.ma.array(np.load(npyfil))
    var[var>1000]=np.ma.masked
    # ISOLATE CCS DOMAIN

    # PLOT TIME SERIES
    for ne in range(var.shape[0]):
        ne_ts = np.mean(np.mean(var[ne,:].squeeze(),axis=2),axis=1)
        if an_val==0:
           ne_ts = [np.mean(ne_ts[12*a:12*(a+1)]) for a in range(nyrs)]
 
        ax[nv].plot(cesm_dates,ne_ts,color='lightgrey')

    em_avg = np.mean(np.mean(np.mean((var),axis=3),axis=2),axis=0).squeeze()

    if an_val == 0: 
       em_avg = [np.mean(em_avg[12*a:12*(a+1)]) for a in range(nyrs)]
  
    ax[nv].plot(cesm_dates,em_avg,color = 'k') 
    ax[nv].set_yticks(ytic[nv])
    plt.setp(ax[nv].get_yticklabels(), fontsize=12)
    nv+=1
ax[0].xaxis_date()
ax[0].set_xlim(dt.datetime(1950,1,1),dt.datetime(2100,12,31))

cesm_ticks = [pltd.date2num(dt.datetime(yr,1,1)) for yr in range(1950,2100+1,25)]
ax[0].set_xticks(cesm_ticks)
plt.setp(ax[2].get_xticklabels(), fontsize=12)
plt.show()
    #print var.shape
    #plt.pcolor(lon,lat,var[0,0,:].squeeze())
    #plt.show()
