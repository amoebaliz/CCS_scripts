
all: regrid_atmos.so

regrid_atmos.so:
	f2py -c -m --fcompiler=intelem --opt='-heap-arrays -fp-model precise -ip -O3 -free -free' regrid_atmos regrid_atmos.f90

