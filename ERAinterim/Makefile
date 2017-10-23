
all: regrid_atmos.so

regrid_atmos.so:
	f2py -c -m --fcompiler=intel --opt='-heap-arrays -fp-model precise -ip -O3 -free -free' regrid_atmos regrid_atmos.f90
#	f2py -c -m --opt='-frepack-arrays -O3 -ffast-math -ffloat-store' regrid_atmos regrid_atmos.f90
#	f2py -c -m  regrid_atmos regrid_atmos.f90
#	f2py -c -m --opt='-frepack-arrays -O3 -ffast-math -ffloat-store -fno-finite-math-only -fmath-errno -ftrapping-math  -frounding-math -ffpe-trap=invalid,zero -ffree-line-length-none' regrid_atmos regrid_atmos.f90
# f2py -c -m --opt='-frepack-arrays' regrid_atmos regrid_atmos.f90
#f2py -c -m --opt='-heap-arrays' regrid_atmos regrid_atmos.f90
#	f2py -c -m --fcompiler=intelem --opt='-heap-arrays -fp-model precise -ip -O3 -free -free' regrid_atmos regrid_atmos.f90

