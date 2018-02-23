
all: regrid_atmos.so

regrid_atmos.so:
	f2py -c -m --fcompiler=intelem --opt='-heap-arrays -fp-model precise -ip -O3 -free -free' regrid_atmos regrid_atmos.f90
# GNU FORTRAN
#	f2py -c -m --opt='-frepack-arrays -O3 -ffast-math -ffloat-store' regrid_atmos regrid_atmos.f90
# ORIGINAL
#	f2py -c -m  regrid_atmos regrid_atmos.f90

