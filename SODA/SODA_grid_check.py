import pyroms
import pyroms_toolbox


soda_grd = '../SODA3.4.1_0.25deg_grid.nc'
srcgrd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(soda_grd, name='SODA3.4.1_0.25deg_CCS',xrange=(1,198),yrange=(1,198))
pyroms_toolbox.BGrid_GFDL.make_remap_grid_file(srcgrd, Bpos='t')

