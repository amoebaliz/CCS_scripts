import os
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4
import pyroms
# -- Adapted from ESMG ROMS Grid generation routine for MOM6 -- #


# Grid dimension

Lm = 382
Mm = 858



# The domain should have right angles at the boundary 
# the LCC and Mercator conserve right angles

lon0 = -142.987301504 ; lat0 = 43.716189653  # upper left
lon1 = -123.001816772 ; lat1 = 17.3303846388  # lower left
lon2 = -109.9 ; lat2 = 22.9  # lower right 
lon3 = -126.7 ; lat3 = 51.2  # upper right

# these are used to create the projection, should contain the corners of the grid.
# See the page for basemap in the matplotlib 

# llcrnrlon = longitude of lower left hand corner of the desired map domain (degrees)   
# llcrnrlat = latitude of lower left hand corner of the desired map domain (degrees)
# urcrnrlon = longitude of upper right hand corner of the desired map domain (degrees)
# urcrnrlat = latitude of upper right hand corner of the desired map domain (degrees)

llcrnrlon = -145
llcrnrlat = 14
urcrnrlon = -109
urcrnrlat = 52

map = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat, resolution = 'l')

map.drawcoastlines()

# corners lon and lat
lonp = np.array([lon0, lon1, lon2, lon3])
latp = np.array([lat0, lat1, lat2, lat3])

# rotation (1 = counter-clockwise, -1 = clockwise)
beta = np.array([1,1,1,1])


# generate the new grid
hgrd = pyroms.grid.Gridgen(lonp,latp, beta,(Mm+3, Lm+3),proj = map)

lonv, latv = map(hgrd.x_vert, hgrd.y_vert, inverse = True)
hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)




plt.show()





# NOTE: attempting to generate grid without PYROMS or PYROMS-TOOLBOX

# PYROMS COMMANDS TO RECODE:

# generate the new grid
#hgrd = pyroms.grid.Gridgen(lonp,latp, beta,(Mm+3, Lm+3),proj = map)
#hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)

# generate the mask
#for verts in map.coastsegs:
#        hgrd.mask_polygon(verts)

# Edit the land mask interactively
#pyroms.grid.edit_mask_mesh(hgrd, proj = map)

# fix minimum depth (here 15 meters)
#hmin = 5
#topo = pyroms_toolbox.change(topo, '<', hmin, hmin)

# ensure that depth is always deeper than hmin
#h = pyroms_toolbox.change(h,'<',hmin,hmin)

# shapiro filter everywhere
#hf = pyroms_toolbox.shapiro_filter.shapiro2(np.log(h),2,1)
#rv = pyroms_toolbox.rvalue(hsharp)
#print 'Max r-value is: ', rv.max()

######## Create the Vertical Grid #################

# vertical coordinate
# defined: 
#theta_b = 2
#theta_s = 7.0
#Tcline = 250
#N = 50
#vgrd = pyroms.vgrid.s_coordinate(h,theta_b, theta_s, Tcline, N, hraw=hraw)

#ROMS grid
#grd_name = 'CCS'
#grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

########### If you need to apply mask changes:
#pyroms.utility.apply_mask_change('mask_change_test.txt',grd)

#########  write grid to netcdf file
#pyroms.grid.write_ROMS_grid(grd, filename = 'new_coarse_grd.nc')


print 'MEEP'


