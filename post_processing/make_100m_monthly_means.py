import numpy as np
import netCDF4 as nc
from os import listdir
from os.path import isfile, join
import glob

# z_r function from Pyroms vgrid
# NOTE: assuming Vtrans==2 even though stated 4 in gridid.txt
def z_r_calc(h, hc, N, s_rho, Cs_r, zeta, Vtrans):
        z_r = np.empty(s_rho.shape + h.shape, 'd')

        for  k in range(N):
             z0 = (hc * s_rho[k] + h * Cs_r[k]) / \
                  (hc + h)
             z_r[k,:] = zeta + (zeta + h * z0)

        return z_r



def z_w_calc(h, hc, Np, s_w, Cs_w, zeta, Vtrans):

        z_w = np.empty((Np) + h.shape, 'd')

        for  k in range(Np):
            z0 = (hc * s_w[k] + h * Cs_w[k]) / \
                          (hc + h)
                    z_w[n,k,:] = zeta[n,:] + (zeta[n,:] + h * z0)

        return z_w

# zslice function from Pyroms
def zslice(var, depth, Cpos='rho', mode='linear'):
    """ 
    zslice = zslice(var, depth, grd)
    optional switch:
      - Cpos='rho', 'u', 'v'  specify the C-grid position where 
			      the variable rely
      - vert=True/False       If True, return the position of 
                              the verticies
      - mode='linear' or 'spline'    specify the type of interpolation
    return a constant-z slice at depth depth from 3D variable var
    lon and lat contain the C-grid position of the slice for plotting.
    If vert=True, lon and lat contain the position of the verticies 
    (to be used with pcolor)
    """
  
    if mode=='linear':
        imode=0
    elif mode=='spline':
        imode=1
    else:
        imode=0
        raise Warning, '%s not supported, defaulting to linear' % mode

    Vtrans = 2
    # N = 50 
    Np = 51
    z_r = z_r_calc(h, hc, 50, s_rho, Cs_r, zeta, 2)
    print z_r.shape
    # compute the depth on Arakawa-C grid position
    
    if Cpos is 'u':
        # average z_r at Arakawa-C u points
        z = 0.5 * (z_r[:,:,:-1] + z_r[:,:,1:])
        mask = mask_u[:]

    elif Cpos is 'v':
        # average z_r at Arakawa-C v points
        z = 0.5 * (z_r[:,:-1,:] + z_r[:,1:,:])
        mask = mask_v[:]

    elif Cpos is 'rho':
        # for temp, salt, rho
        z = z_r[:]
        mask = mask_rho[:]

    else:
        raise Warning, '%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos

    assert len(z.shape) == 3, 'z must be 3D'
    assert len(var.shape) == 3, 'var must be 3D'
    assert z.shape == var.shape, 'data and prop must be the same size'

    depth = -abs(depth)
    depth = depth * np.ones(z.shape[1:])
    
    zslice = _iso.zslice(z, var, imode)

    # mask land
    zslice = np.ma.masked_where(mask == 0, zslice)
    # mask region with shalower depth than requisted depth
    zslice = np.ma.masked_where(zslice == 1e20, zslice)

    return zslice

# ROMS Grid information
grdfile = '/Users/elizabethdrenkard/ANALYSES/CCS/Inputs/Grid/CCS_grd_high_res_bathy_jerlov.nc'
grdfid = nc.Dataset(grdfile)
mask = grdfid.variables['mask_rho'][:]
lat_rho = grdfid.variables['lat_rho'][:]
lon_rho = grdfid.variables['lon_rho'][:]
lat_u = grdfid.variables['lat_u'][:]
lon_u = grdfid.variables['lon_u'][:]
lat_v = grdfid.variables['lat_v'][:]
lon_v = grdfid.variables['lon_v'][:]
lat_vert = grdfid.variables['lat_vert'][:]
lon_vert = grdfid.variables['lon_vert'][:]
mask_rho = grdfid.variables['mask_rho'][:]
mask_u = grdfid.variables['mask_u'][:]
mask_v = grdfid.variables['mask_v'][:]
h = grdfid.variables['h'][:]
hc = grdfid.variables['hc'][:]
s_rho = grdfid.variables['s_rho'][:]
Cs_r = grdfid.variables['Cs_r'][:]

# Depth of interest meters
depth=100


for yr in range(3,12+1):
    #print yr

    for mn in range(12):
        #print mn+1
        mypath = '/Volumes/Abalone/CCS/his2/' + str(yr).zfill(4) + \
                 '/CCS-LD.HCo02Y_avg_' + str(yr).zfill(4) + \
                 '-' + str(mn+1).zfill(2) + '*T12:00:00.nc'
           
        # list all files of given year-month make sure they exist. 
        # NOTE: day order in month will not be preserved BUT, this 
        #       should be unnecessary if commputing monthly means
        n=0
        for f in glob.glob(mypath):
            print f 
            fid = nc.Dataset(f)
            zeta = fid.variables['zeta'][:].squeeze() 
            temp = fid.variables['temp'][:].squeeze()
            u = fid.variables['u'][:].squeeze()
            v = fid.variables['v'][:].squeeze()
            
            # Initalize storage fields
            if n==0:
               temp_stor = np.zeros(lat_rho.shape)
               u_stor = np.zeros(lat_u.shape)
               v_stor =  np.zeros(lat_v.shape)

            # UPDATE STORAGE FIELDS WITH CALCULATED VALUE AT DESIRED DEPTH
            # NOTE: Using default linear interpolation (as opposed to a spine)
            temp_stor+= zslice(temp, depth, Cpos='rho')
            u_stor+= zslice(u, depth, Cpos='u')
            v_stor+= zslice(v, depth, Cpos='v')
 
            n+=1

        # Generate monthly mean fields by dividing by n
        temp_depth = temp_stor/n
        u_depth = u_stor/n
        v_depth = v_stor/n

        # Save as monthly netCDF file



