import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import ESMF
import shapefile as shp 

mon_list =['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'] 
filname = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/collaborators/Farrah_shp/caltrawl_GCS'     
r = shp.Reader(filname)

w = shp.Writer(r.shapeType)
w.fields = [r.fields[n] for n in range(2)]
for nm in range(12):
    w.fields.append([mon_list[nm]] + ['N', 10, 5])

w.records = []
months = list(np.arange(12)+1)

for nr in range(len(r.records())):
    w.records.append(r.records()[nr][:1] + months)    
w._shapes.extend(r.shapes()) 

w.save('test')



