import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import shapefile as shp


# ROMS Grid information
grdfile = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
fid = nc.Dataset(grdfile)
mask_rho = fid.variables['mask_rho'][:]
rlat = fid.variables['lat_rho'][:]
rlon = fid.variables['lon_rho'][:]
plat = fid.variables['lat_psi'][:]
plon = fid.variables['lon_psi'][:]
dep = fid.variables['h'][:]

# WEB: https://pypi.org/project/pyshp/

shp_fil = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/collaborators/Farrah_shp/caltrawl_GCS'

sf = shp.Reader(shp_fil)
nshapes = len(sf.shapes())

# elements in the records
sf.fields

# records corresponding to the fields
sf.records()

m_offset = 0.05
mask_val = 0
map_order = 30

plot_num = 2 

lat_bnds = np.array(((32.35,34.49),  \
                     (34.51,37.99),\
                     (38.01,42)))

lon_bnds = np.array(((-122.05,-117.065),  \
                     (-124.4,-120.455),\
                     (-125.7,-122.8)))

fig_sz = np.array(((10,4),(8,7),(6,8)))

fig, ax = plt.subplots(figsize=(2*abs(np.diff(lon_bnds[plot_num,:])),2*abs(np.diff(lat_bnds[plot_num,:]))))

m = Basemap(llcrnrlat=np.min(plat)-m_offset,urcrnrlat = np.max(plat)+m_offset,llcrnrlon=np.min(plon)-m_offset,urcrnrlon=np.max(plon)+m_offset, resolution='i', ax=ax)

for ns in range(nshapes):
    s = sf.shape(ns)
    bbox = [coord for coord in s.bbox]
    # EXTRACT CORNER COORDINATES
    shp_lons = np.array((bbox[0], bbox[2]))
    shp_lats = np.array((bbox[3], bbox[1]))

    # IF INSIDE MAP SUBDOMAIN
    if ((lon_bnds[plot_num,0]<shp_lons[0]<lon_bnds[plot_num,1])  or  \
        (lon_bnds[plot_num,0]<shp_lons[1]<lon_bnds[plot_num,1])) and \
       ((lat_bnds[plot_num,0]<shp_lats[0]<lat_bnds[plot_num,1])  or  \
        (lat_bnds[plot_num,0]<shp_lats[1]<lat_bnds[plot_num,1])):
 
        x = [bbox[0], bbox[2], bbox[2], bbox[0], bbox[0]]
        y = [bbox[3], bbox[3], bbox[1], bbox[1], bbox[3]]
        m.plot(x,y)

ax.set_ylim(lat_bnds[plot_num,0]-m_offset,lat_bnds[plot_num,1]+m_offset)
ax.set_xlim(lon_bnds[plot_num,0],lon_bnds[plot_num,1])

tx =plt.text(np.max(lon_bnds[plot_num,:])-.5,\
             np.max(lat_bnds[plot_num,:])-.25,\
             'MAY', fontsize=20,zorder=map_order+5)

m.drawmeridians(np.arange(-125,-119+1,2), labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)
m.drawparallels(np.arange(34,42+1,2), labels=[0,0,0,0], fmt='%d', fontsize=18,zorder=map_order+5)

plt.show()

ns=0
for shape in sf.shapeRecords():
    s = sf.shape(ns)
    # Round coordinates to 3 decimal places
#    print ['%.3f' % coord for coord in s.bbox]    
#    print [i[0] for i in shape.shape.points[:]]
#    print [i[1] for i in shape.shape.points[:]]
#    print
#    
    ns+=1 






# LOOP OVER EACH SHAPE
for ns in range(nshapes):
    s = sf.shape(ns)
    bbox = [coord for coord in s.bbox]
    #x = [bbox[0], bbox[2], bbox[2], bbox[0], bbox[0]]
    #y = [bbox[3], bbox[3], bbox[1], bbox[1], bbox[3]]
    #plt.plot(x,y)

    c_lon = np.mean((bbox[0], bbox[2]))
    c_lat = np.mean((bbox[3], bbox[1]))

    x = [2*bbox[0]-c_lon, c_lon, 2*bbox[2]-c_lon]
    y = [2*bbox[1]-c_lat, c_lat, 2*bbox[3]-c_lat]

    Xn, Yn = np.meshgrid(x, y)
 
    #for nx in x:
    #    for ny in y: 
    #        plt.plot(nx,ny,'ob')

    #plt.plot(c_lon,c_lat,'ok')
    #plt.show()
# WRITE SHAPE FILES FOR... SST

# ADD RECORDS FOR ... CHANGE IN SST
