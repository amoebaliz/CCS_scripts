import numpy as np
import netCDF4 as nc
import csv
import ESMF
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import datetime as dt
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection

def polygon_patch(mapid,axs):
    mapid.drawcoastlines(linewidth=0,zorder=map_order+6)
    mapid.drawmapboundary(fill_color=[1,1,1])
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor=(1,1,1), closed=False)
    axs.add_collection(lc)


def ccs_region(lat,lon):
    # CCS ROMS GRID FOR ISOLATING CCS
    row_min = 300
    row_max = -370
    col_min = 150
    fid = nc.Dataset('/glade/work/edrenkar/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc')
    r_lat = fid.variables['lat_rho'][row_min:row_max,col_min:]
    r_lon = fid.variables['lon_rho'][row_min:row_max,col_min:]
    sourcegrid = ESMF.Grid(np.array(r_lon.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)
    # CESM GRID
    lon[lon>180]=lon[lon>180]-360
    destgrid = ESMF.Grid(np.array(lon.shape), staggerloc = ESMF.StaggerLoc.CENTER, coord_sys = ESMF.CoordSys.SPH_DEG)

    ## POINTERS
    source_lon = sourcegrid.get_coords(0)
    source_lat = sourcegrid.get_coords(1)
    dest_lon = destgrid.get_coords(0)
    dest_lat = destgrid.get_coords(1)
    ## FILLS
    source_lon[...] = r_lon
    source_lat[...] = r_lat
    dest_lon[...] = lon
    dest_lat[...] = lat

    sourcefield = ESMF.Field(sourcegrid, name = 'CCS ROMS')
    destfield = ESMF.Field(destgrid, name = 'CESM Sub')

    sourcefield.data[...] = fid.variables['mask_rho'][row_min:row_max,col_min:].squeeze()
    regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.BILINEAR,
             unmapped_action = ESMF.UnmappedAction.IGNORE)
    destfield = regrid(sourcefield, destfield)

    return_val = destfield.data

    # CONVERTING TO 1/0 boolean for mask_where: TRUE (1) where you want to mask
    return_val[return_val<.5]=0    # FALSE (initially land and ocean outside ROMS domain) 
    return_val[return_val>.5]=1    # TRUE  (initially water inside ROMS domain)

    return_val = -1*return_val + 1 # REVERSE MASK VALUES FOR where_mask 
    return_val2 = return_val.astype(np.bool) # MAKE BOOLEAN

    return return_val2


# OBJECTIVE: Plot Surf Temp, O2, PH time-series
dir = '/glade/work/edrenkar/'
grid_fil = dir + 'LENS_CCS_lat_lon.npy'
station_fil = 'CalCOFIStaPosNDepth113.csv'

latlon = np.load(grid_fil)
lat = latlon[0,:].squeeze()
lon = latlon[1,:].squeeze()

CCS_MASK = ccs_region(lat,lon)

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
fig, ax = plt.subplots(2,2, figsize=(9,7))
ytic = [range(15,23+1,3), np.arange(7.7,8.1+.1,.2), range(235,265+1,15)]        
ytit = ['Temperature ($^o$C)', 'pH', 'Oxygen (mmol/m$^3$)']
ystr = ['b','c','d']
vary = [21,8.13,262.5]
cesm_ticks = [pltd.date2num(dt.datetime(yr,1,1)) for yr in range(1950,2100+1,25)]

nv=0
em_str = np.zeros((nyrs,3))
for VAR in ("SST", "PH", "O2"):
    npyfil = dir + 'LENS_1950-2100_' + VAR + '.npy'
    var = np.ma.array(np.load(npyfil))
    var[var>1000]=np.ma.masked
    # ISOLATE CCS DOMAIN
    var[:,:,CCS_MASK] = np.ma.masked
    if nv==0:
       var_dom = np.ma.where(var>0,3,1)
       var_dom = np.ma.masked_array(var_dom, var_dom == 1)
       # GET GEO Vertices
       lat_vert = (lat[:-1,:-1] + lat[:-1,1:] + lat[1:,:-1] + lat[1:,1:])/4.
       lon_vert = (lon[:-1,:-1] + lon[:-1,1:] + lon[1:,:-1] + lon[1:,1:])/4.
       # Draw map domain
       map_order=30
       m_lat = [31,37]; m_lon = [-126,-118]; m_off=1.5
       m = Basemap(llcrnrlat=m_lat[0]-m_off,urcrnrlat=m_lat[1]+m_off,\
           llcrnrlon=m_lon[0]-m_off,urcrnrlon=m_lon[1]+m_off, resolution='i',ax=ax[0,0])
       P = m.pcolormesh(lon_vert,lat_vert,np.mean(np.mean(var_dom,axis=1),axis=0)[1:-1,1:-1],vmin=0,vmax=10,cmap='Greys')
       P2 = m.pcolormesh(lon_vert,lat_vert,np.ma.mean(np.ma.mean(var_dom,axis=1),axis=0)[1:-1,1:-1],facecolor='none',edgecolors=(1.0, 1.0, 1.0, 0.3),linewidth=0.005)
       P2._is_stroked = False
       # NOTE: I don't understand why 'facecolor=none' or particular edgecolors values are mandatory
       polygon_patch(m,ax[0,0])
       # plot CalCOFI stations
       with open(station_fil) as f:
            csv_reader = csv.reader(f)
            next(csv_reader, None) # skip first row of headers
            for row in csv_reader:
                m.plot(-1*float(row[4]),float(row[3]),'ko',ms=2)

       m.drawmeridians(m_lon, labels=[0,0,1,0], fmt='%d', fontsize=12)
       m.drawparallels(m_lat, labels=[1,0,0,0], fmt='%d', fontsize=12)
       ax[0,0].text(-127, 37.5, 'a', fontsize=14) 
       #plt.show()

    # INDICES FOR PLOTS
    rv = (nv+1)/2
    cv = abs(-1*nv+1) 
    # PLOT TIME SERIES
    for ne in range(var.shape[0]):
        ne_ts = np.ma.mean(np.ma.mean(var[ne,:].squeeze(),axis=2),axis=1)
        if an_val==0:
           ne_ts = [np.ma.mean(ne_ts[12*a:12*(a+1)]) for a in range(nyrs)]
 
        ax[rv,cv].plot(cesm_dates,ne_ts,color='lightgrey')

    em_avg = np.ma.mean(np.ma.mean(np.ma.mean((var),axis=3),axis=2),axis=0).squeeze()
    if an_val == 0: 
       em_avg = [np.ma.mean(em_avg[12*a:12*(a+1)]) for a in range(nyrs)]

    em_str[:,nv] = em_avg 
    ax[rv,cv].plot(cesm_dates,em_avg,color = 'k')
    
    if nv==1:
       ax[rv,cv].set_ylim(7.7,8.17) 
    ax[rv,cv].set_yticks(ytic[nv])
    ax[rv,cv].set_ylabel(ytit[nv],fontsize=12)
    ax[rv,cv].text(dt.datetime(1953,1,1), vary[nv], ystr[nv], fontsize=14)
    plt.setp(ax[rv,cv].get_yticklabels(), fontsize=12)
    ax[rv,cv].xaxis_date()
    ax[rv,cv].set_xlim(dt.datetime(1950,1,1),dt.datetime(2100,12,31))
    ax[rv,cv].set_xticks(cesm_ticks)
    plt.setp(ax[rv,cv].get_xticklabels(), fontsize=11)
    nv+=1
plt.tight_layout()
plt.savefig('Annual_mean_calcofi_domain')

# sfname = 'annual_mean_SST_O2_pH'
# SAVE TO FILE
# with open(sfname, 'w') as f:
#     writer = csv.writer(f, delimiter=' ', lineterminator='\n')
#     for n in range(nyrs):
#         row = [str(1950+n)] + [str(em_str[n,0])] + [str(em_str[n,1])] + [str(em_str[n,2])] 
#         writer.writerow(row)

plt.show()
