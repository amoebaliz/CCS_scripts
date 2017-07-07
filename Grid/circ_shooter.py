from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
 
import numpy as np
 
### PARAMETERS FOR MATPLOTLIB :
import matplotlib as mpl
mpl.rcParams['font.size'] = 10.
mpl.rcParams['font.family'] = 'Comic Sans MS'
mpl.rcParams['axes.labelsize'] = 8.
mpl.rcParams['xtick.labelsize'] = 6.
mpl.rcParams['ytick.labelsize'] = 6.
 
def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.
 
    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")
 
    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)
 
    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):
 
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu
 
    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
 
    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
 
    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi
 
    return (glon2, glat2, baz)
 
def great(m, startlon, startlat, azimuth,*args, **kwargs):
    glon1 = startlon
    glat1 = startlat
    glon2 = glon1
    glat2 = glat1
 
    step = 50
 
    glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
    if azimuth-180 >= 0:
        while glon2 <= startlon:
            m.drawgreatcircle(glon1, glat1, glon2, glat2,del_s=50,**kwargs)
            azimuth = baz + 180.
            glat1, glon1 = (glat2, glon2)
 
            glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
    elif azimuth-180 < 0:
        while glon2 >= startlon:
            m.drawgreatcircle(glon1, glat1, glon2, glat2,del_s=50,**kwargs)
            azimuth = baz + 180.
            glat1, glon1 = (glat2, glon2)
 
            glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
 
fig = plt.figure(figsize=(11.7,8.3))
 
#Custom adjust of the subplots
 
plt.subplots_adjust(left=0.05,right=0.95,top=0.90,bottom=0.05,wspace=0.15,hspace=0.05)
 
ax = plt.subplot(111)
 
#Let's create a basemap of the world
llcrnrlon = -150
llcrnrlat = 14
urcrnrlon = -109
urcrnrlat = 52
 
m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat, resolution = 'l')
m.drawcountries()
m.drawcoastlines()
 
glon1 = -109.9
glat1 = 22.9
 
azimuth = 339
maxdist = 3450

glon2, glat2, baz = shoot(glon1, glat1, azimuth, maxdist)

m.drawgreatcircle(glon1, glat1, glon2, glat2,del_s=50,color='red', lw=1.)
print glon2, glat2

#gc_pts = np.array(m.gcpoints(glon1, glat1, glon2, glat2, 3))
#lonpt, latpt = m(gc_pts[0,1],gc_pts[1,1],inverse=True)
#print lonpt,latpt
#m.plot(gc_pts[0,:],gc_pts[1,:],'ob')
#m.plot(lonpt,latpt,'og')
#m.plot(gc_pts[0,1],gc_pts[1,1],'og')

azimuth = 340-90.
maxdist = 1500

#glon_C, glat_C, baz = shoot(lonpt, latpt, azimuth, maxdist)
#m.drawgreatcircle(lonpt, latpt, glon_C, glat_C,del_s=50,color='blue', lw=4.)

glon_L, glat_L, baz = shoot(glon1, glat1, azimuth, maxdist)
glon_U, glat_U, baz = shoot(glon2, glat2, azimuth, maxdist)

print 'lower right', glon1, glat1 
print 'upper right', glon2, glat2
print 'lower left', glon_L, glat_L
print 'upper left' ,glon_U, glat_U

m.drawgreatcircle(glon2, glat2, glon_U, glat_U,del_s=50,color='green', lw=1.)
m.drawgreatcircle(glon1, glat1, glon_L, glat_L,del_s=50,color='orange', lw=1.) 
m.drawgreatcircle(glon_U, glat_U, glon_L, glat_L,del_s=50,color='blue', lw=1.) 

#m.drawgreatcircle(glon1, glat1, glon_C, glat_C,del_s=50,color='green', lw=4.)
#m.drawgreatcircle(glon2, glat2, glon_C, glat_C,del_s=50,color='orange', lw=4.)
 
#startlon = 4.360515
#startlat = 50.79747
#azimuth = 165.
#great(m,startlon, startlat, azimuth, color='orange', lw=2.0)
#great(m,startlon, startlat, azimuth+180., color='orange', lw=2.0)
 
#plt.savefig('tutorial08.png',dpi=300)
 
plt.show()

