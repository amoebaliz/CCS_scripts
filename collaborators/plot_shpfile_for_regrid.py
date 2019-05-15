import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import shapefile as shp

# WEB: https://pypi.org/project/pyshp/

shp_fil = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/colaborators/Farrah_shp/caltrawl_GCS'

sf = shp.Reader(shp_fil)
nshapes = len(sf.shapes())

# elements in the records
sf.fields

# records corresponding to the fields
sf.records()

plt.figure()
for ns in range(nshapes):
    s = sf.shape(ns)
    bbox = [coord for coord in s.bbox] 
    x = [bbox[0], bbox[2], bbox[2], bbox[0], bbox[0]]
    y = [bbox[3], bbox[3], bbox[1], bbox[1], bbox[3]]
    plt.plot(x,y)
#plt.show()

ns=0
for shape in sf.shapeRecords():
    s = sf.shape(ns)
    # Round coordinates to 3 decimal places
    print ['%.3f' % coord for coord in s.bbox]    
    print [i[0] for i in shape.shape.points[:]]
    print [i[1] for i in shape.shape.points[:]]
    print
#    
    ns+=1 


plt.figure()
for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
#    print x, y
    plt.plot(x,y)
plt.show()


plt.figure()





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
