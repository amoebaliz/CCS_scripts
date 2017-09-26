import netCDF4 as nc
import numpy as np
import calendar
import datetime as dt
import humidity_toolbox
import os 

class ERAinterim_processing():
	''' A Class to perform operation on ERAinterim files 
	such as decumulation, means, unit conversions,... '''

	def __init__(self,dict_input):
		# constants
		self.rho_w = 1000.
		self.nsec_per_day = 86400.
		self.reftime = dt.datetime(1900,1,1,0,0)
		self.spval = 1.0e+15
		self.dataset = 'ERAinterim'
		# read inputs
		self.dict_input = dict_input
		for key in dict_input:
                	exec('self.' + key + '=dict_input[key]')
		# set time related stuff
		#if calendar.isleap(self.year):
		#	self.ndays = 366
		#else:
		#	self.ndays = 365
		# forcing set freq
		self.nframes_per_day = 8

		self.nframes = self.ndays * self.nframes_per_day
		
		print self.year, self.month, 'has', self.nframes, 'frames'

		return None

	def __call__(self):
		if self.dict_input.has_key('file_precip'):
			print 'Processing precip file...'
			self.process_precip_to_daily()
		if self.dict_input.has_key('file_radlw'):
			print 'Processing longwave file...'
			self.process_radlw_to_daily()
		if self.dict_input.has_key('file_radsw'):
			print 'Processing shortwave file...'
			self.process_radsw_to_daily()
		if self.dict_input.has_key('file_d2'):
			print 'Create specific humidity file...'
			self.create_q2_file()
		if self.dict_input.has_key('file_t2'):
			print 'Rewrite t2 file...'
			self.process_t2_file()
		if self.dict_input.has_key('file_msl'):
			print 'Rewrite msl file...'
			self.process_msl_file()
		if self.dict_input.has_key('file_u10'):
			print 'Rewrite u10 file...'
			self.process_u10_file()
		if self.dict_input.has_key('file_v10'):
			print 'Rewrite v10 file...'
			self.process_v10_file()
		return None

	#------------------ Meta functions ------------------------------------------
		
	def process_precip_to_daily(self):
		''' processing the precipitation file :
		sum all cumulated fields and divide by a day'''
		precip_out = np.empty((self.ndays,self.ny,self.nx))
		time = np.empty((self.ndays))
		# open file
		fid_precip = self._opennc(self.file_precip)
		# read coordinates and time
		lon = self._readnc(fid_precip,'lon')
		lat = self._readnc(fid_precip,'lat')
		# run the summation
		for kt in np.arange(0,self.ndays):
			tmp = np.zeros((self.ny,self.nx))
			nvalues = self.nframes_per_day / self.ncumul # number of values to read
			for kc in np.arange(nvalues):
				# frames to read
				# ERAinterim (nframes_per_day = 8) : for kt = 0 (day 1) we read frames 4 and 8
				kframe = (kt * self.nframes_per_day) + (kc+1) * self.ncumul
				tmp = tmp + self._readnc_oneframe(fid_precip,'TP',kframe-1) # C indexing hence -1

			precip_out[kt,:,:] = tmp.copy()
			this_day = dt.datetime(self.year,int(self.month),1,12,0) + dt.timedelta(days=int(kt))
			time[kt] = (this_day - self.reftime).days + (this_day - self.reftime).seconds / 86400.
 
		# close file
		self._closenc(fid_precip)
		# units conversion (from m to m/s then to kg/m2/s)
		cumul_time = self.nsec_per_day 
		precip_out = precip_out * self.rho_w / cumul_time
		# remove negative values
		precip_out[np.where(precip_out < 0.)] = 0.
                 
                # create monthly mean
	        precip_out = np.mean(precip_out,axis=0,keepdims=True) 
                time = np.mean(time)
                
                # write file
		my_dict = {'varname':'rain','time_dim':'time','time_var':'time','long name':'Total Precipitation',\
			'units':'kg.m-2.s-1','fileout':self.output_dir + 'precip_' + self.dataset + '_' + str(self.year) + '_' + str(self.month) + '_monthly_ROMS.nc'}
		self._write_ncfile(lon,lat[::-1],time,precip_out[:,::-1,:],my_dict)
		precip_out = None
		return None

        def process_radlw_to_daily(self):
                ''' processing the longwave radiation file :
                sum all cumulated fields and divide by a day'''
                radlw_out = np.empty((self.ndays,self.ny,self.nx))
                time = np.empty((self.ndays))
                # open file
                fid_radlw = self._opennc(self.file_radlw)
                # read coordinates and time
                lon = self._readnc(fid_radlw,'lon')
                lat = self._readnc(fid_radlw,'lat')
                # run the decumulation
                for kt in np.arange(0,self.ndays):
                        tmp = np.zeros((self.ny,self.nx))
                        nvalues = self.nframes_per_day / self.ncumul # number of values to read
                        for kc in np.arange(nvalues):
                                # frames to read
                                # ERAinterim (nframes_per_day = 8) : for kt = 0 (day 1) we read frames 4 and 8
                                kframe = (kt * self.nframes_per_day) + (kc+1) * self.ncumul
                                tmp = tmp + self._readnc_oneframe(fid_radlw,'STRD',kframe-1) # C indexing hence -1
                        radlw_out[kt,:,:] = tmp.copy()
                        this_day = dt.datetime(self.year,int(self.month),1,12,0) + dt.timedelta(days=int(kt))
                        time[kt] = (this_day - self.reftime).days + (this_day - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_radlw)
                # units conversion (from W.m-2.s to W.m-2)
                cumul_time = self.nsec_per_day
                radlw_out = radlw_out / cumul_time

                # create monthly mean
                radlw_out = np.mean(radlw_out,axis=0,keepdims=True)
                time = np.mean(time)

                # write file
	        my_dict = {'varname':'lwrad_down','time_dim':'time','time_var':'time','long name':'Downwelling longwave radiation',\
	                'units':'W.m-2','fileout':self.output_dir + 'radlw_' + self.dataset + '_' + str(self.year) + '_' + str(self.month) + '_monthly_ROMS.nc'}
       		self._write_ncfile(lon,lat[::-1],time,radlw_out[:,::-1,:],my_dict)
		radlw_out = None
                return None

        def process_radsw_to_daily(self):
                ''' processing the shortwave radiation file :
                sum all 3h fields and divide by a day'''
                radsw_out = np.empty((self.ndays,self.ny,self.nx))
                time = np.empty((self.ndays))
                # open file
                fid_radsw = self._opennc(self.file_radsw)
                # read coordinates and time
                lon = self._readnc(fid_radsw,'lon')
                lat = self._readnc(fid_radsw,'lat')
                # run the summation
                for kt in np.arange(0,self.ndays):
                        tmp = np.zeros((self.ny,self.nx))
                        nvalues = self.nframes_per_day / self.ncumul # number of values to read
                        for kc in np.arange(nvalues):
                                # frames to read
                                # ERAinterim (nframes_per_day = 8) : for kt = 0 (day 1) we read frames 4 and 8
                                kframe = (kt * self.nframes_per_day) + (kc+1) * self.ncumul
                                tmp = tmp + self._readnc_oneframe(fid_radsw,'SSRD',kframe-1) # C indexing hence -1
                        radsw_out[kt,:,:] = tmp.copy()
                        this_day = dt.datetime(self.year,int(self.month),1,12,0) + dt.timedelta(days=int(kt))
                        time[kt] = (this_day - self.reftime).days + (this_day - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_radsw)
                # units conversion (from W.m-2.s to W.m-2)
                cumul_time = self.nsec_per_day
                radsw_out = radsw_out / cumul_time

                # create monthly mean
                radsw_out = np.mean(radsw_out,axis=0,keepdims=True)
                time = np.mean(time)
                
                # write file
	        my_dict = {'varname':'swrad','time_dim':'time','time_var':'time','long name':'Shortwave radiation',\
	                'units':'W.m-2','fileout':self.output_dir + 'radsw_' + self.dataset + '_' + str(self.year) + '_' + str(self.month) + '_monthly_ROMS.nc'}
	        self._write_ncfile(lon,lat[::-1],time,radsw_out[:,::-1,:],my_dict)
		radsw_out = None
                return None

	def create_q2_file(self):
		''' Create specific humidity file '''
		q2_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_d2 = self._opennc(self.file_d2)
                fid_msl = self._opennc(self.file_msl)
                # read coordinates and time
                lon = self._readnc(fid_d2,'lon')
                lat = self._readnc(fid_d2,'lat')
		# run the computation
		for kt in np.arange(0,self.nframes):
			d2 = self._readnc_oneframe(fid_d2,'D2M',kt)
			msl = self._readnc_oneframe(fid_msl,'MSL',kt)
			q2_tmp = humidity_toolbox.q2_from_d2_and_msl(d2,msl)
			q2_out[kt,:,:] = q2_tmp.copy()
			this_time = dt.datetime(self.year,int(self.month),1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_d2)
                self._closenc(fid_msl)

                # create monthly mean
                q2_out = np.mean(q2_out,axis=0,keepdims=True) 
                time = np.mean(time)
                # write file
	        my_dict = {'varname':'Qair','time_dim':'time','time_var':'time','long name':'Specific humidity at 2m',\
	                'units':'kg/kg','fileout':self.output_dir + 'q2_' + self.dataset + '_' + str(self.year) + '_' + str(self.month) + '_monthly_ROMS.nc'}
	        self._write_ncfile(lon,lat[::-1],time,q2_out[:,::-1,:],my_dict)
		q2_out = None
                return None

	def process_t2_file(self):
		''' Rewrite temperature file according to model's needs '''
		t2_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_t2 = self._opennc(self.file_t2)
                # read coordinates and time
                lon = self._readnc(fid_t2,'lon')
                lat = self._readnc(fid_t2,'lat')
		# run the computation
		for kt in np.arange(0,self.nframes):
			t2_out[kt,:,:] = self._readnc_oneframe(fid_t2,'T2M',kt) - 273.15
			this_time = dt.datetime(self.year,int(self.month),1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_t2)

                # create monthly mean
                t2_out = np.mean(t2_out,axis=0,keepdims=True) 
                time = np.mean(time)

                # write file
	        my_dict = {'varname':'Tair','time_dim':'time','time_var':'time','long name':'Air Temperature at 2m',\
	                'units':'degC','fileout':self.output_dir + 't2_' + self.dataset + '_' + str(self.year) + '_' + str(self.month) + '_monthly_ROMS.nc'}
	        self._write_ncfile(lon,lat[::-1],time,t2_out[:,::-1,:],my_dict)
		t2_out = None
                return None

	def process_msl_file(self):
		''' Rewrite pressure file according to model's needs '''
		msl_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_msl = self._opennc(self.file_msl)
                # read coordinates and time
                lon = self._readnc(fid_msl,'lon')
                lat = self._readnc(fid_msl,'lat')
		# run the computation
		for kt in np.arange(0,self.nframes):
			msl_out[kt,:,:] = self._readnc_oneframe(fid_msl,'MSL',kt)
			this_time = dt.datetime(self.year,int(self.month),1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_msl)

                # create monthly mean
                msl_out = np.mean(msl_out,axis=0,keepdims=True) 
                time = np.mean(time)

                # write file
	        my_dict = {'varname':'Pair','time_dim':'time','time_var':'time','long name':'Mean sea-level pressure',\
	                'units':'Pa','fileout':self.output_dir + 'msl_' + self.dataset + '_' + str(self.year) + '_' + str(self.month) + '_monthly_ROMS.nc'}
	        self._write_ncfile(lon,lat[::-1],time,msl_out[:,::-1,:],my_dict)
		msl_out = None
		return None

	def process_u10_file(self):
		''' Rewrite zonal wind file according to model's needs '''
		u10_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_u10 = self._opennc(self.file_u10)
                # read coordinates and time
                lon = self._readnc(fid_u10,'lon')
                lat = self._readnc(fid_u10,'lat')
		# run the computation
		for kt in np.arange(0,self.nframes):
			u10_out[kt,:,:] = self._readnc_oneframe(fid_u10,'U10M',kt)
			this_time = dt.datetime(self.year,int(self.month),1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_u10)

                # create monthly mean
                u10_out = np.mean(u10_out,axis=0,keepdims=True) 
                time = np.mean(time)

                # write file
	        my_dict = {'varname':'Uwind','time_dim':'time','time_var':'time','long name':'Zonal wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'u10_' + self.dataset + '_' + str(self.year) + '_' + str(self.month) + '_monthly_ROMS.nc'}
	        self._write_ncfile(lon,lat[::-1],time,u10_out[:,::-1,:],my_dict)
		u10_out = None
		return None

	def process_v10_file(self):
		''' Rewrite meridional wind file according to model's needs '''
		v10_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_v10 = self._opennc(self.file_v10)
                # read coordinates and time
                lon = self._readnc(fid_v10,'lon')
                lat = self._readnc(fid_v10,'lat')
		# run the computation
		for kt in np.arange(0,self.nframes):
			v10_out[kt,:,:] = self._readnc_oneframe(fid_v10,'V10M',kt)
			this_time = dt.datetime(self.year,int(self.month),1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_v10)

                # create monthly mean
                v10_out = np.mean(v10_out,axis=0,keepdims=True) 
                time = np.mean(time)

                # write file
	        my_dict = {'varname':'Vwind','time_dim':'time','time_var':'time','long name':'Meridional wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'v10_' + self.dataset + '_' + str(self.year) + '_' + str(self.month) + '_monthly_ROMS.nc'}
	        self._write_ncfile(lon,lat[::-1],time,v10_out[:,::-1,:],my_dict)
		v10_out = None
		return None

	#------------------ NetCDF functions ------------------------------------------
	def _opennc(self,myfile):
		fid = nc.Dataset(myfile,'r')
		return fid

	def _closenc(self,fid):
		fid.close()
		return None

        def _readnc_frames(self,fid,myvar,myframe_start,myframe_end):
                ''' read data from netcdf '''
                out = fid.variables[myvar][myframe_start:myframe_end,:,:].squeeze()
                return out

        def _readnc_oneframe(self,fid,myvar,myframe):
                ''' read data from netcdf '''
                out = fid.variables[myvar][myframe,:,:].squeeze()
                return out

        def _readnc(self,fid,myvar):
                ''' read data from netcdf '''
                out = fid.variables[myvar][:].squeeze()
                return out

	def _write_ncfile(self,lon_array,lat_array,time,var,dict_wrt):
        	fid = nc.Dataset(dict_wrt['fileout'], 'w', format='NETCDF3_CLASSIC')
	        fid.description = 'ERAinterim post-processing based on code by raphael.dussin@gmail.com'
	        # dimensions
	        fid.createDimension('lat', lat_array.shape[0])
	        fid.createDimension('lon', lon_array.shape[0])
	        fid.createDimension(dict_wrt['time_dim'], None)
	        # variables
	        latitudes  = fid.createVariable('lat', 'f8', ('lat',))
	        longitudes = fid.createVariable('lon', 'f8', ('lon',))
	        times      = fid.createVariable(dict_wrt['time_var'], 'f8', (dict_wrt['time_dim'],))
	        variable   = fid.createVariable(dict_wrt['varname'], 'f4', (dict_wrt['time_dim'],'lat','lon',),fill_value=self.spval)
	
		# attributes
		longitudes.units = "degrees_east" 
		longitudes.valid_min = lon_array.min()
		longitudes.valid_max = lon_array.max()
		longitudes.long_name = "longitude" 

		latitudes.units = "degrees_north" 
		latitudes.valid_min = lat_array.min()
		latitudes.valid_max = lat_array.max()
		latitudes.long_name = "latitude" 

		times.units = "days since " + self.reftime.isoformat()
		#times.valid_min = time.min()
                #times.valid_max = time.max()
                times.calendar = "NO LEAP"

		variable.long_name = dict_wrt['long name']
		variable.units = dict_wrt['units']
		variable.coordinates = "lon lat" 
		variable.time = dict_wrt['time_var']
		variable.missing_value = self.spval
		#variable.valid_range = var.min() , var.max()

	        # data
	        latitudes[:]    = lat_array
	        longitudes[:]   = lon_array
	        times[:]        = time
	        variable[:,:,:] = var

	        # close
	        fid.close()
	        return None

