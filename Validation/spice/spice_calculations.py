# SPICE CALCULATION
# adapted from algorithm developed by 
# P. Flament

# http://www.satlab.hawaii.edu/spice/spice.m


import numpy as np

# Generate coefficent matrix, B
spc_fil = '/Users/elizabethdrenkard/TOOLS/CCS_scripts/Validation/spice/spice.txt'


B = np.empty([5,6])
with open(spc_fil) as sfil:
     next(sfil)
     for line in sfil:
         B[int(line.split()[0]),int(line.split()[1])]=float(line.split()[2])*10**int(line.split()[3][1:-1])

# READ IN TEMP and SALT
t = 
s = 

# DIMENSIONALIZE STORAGE VARIABLES
rw, cl = t.shape          # dimensions of temp field
sp = np.zeros([rw,cl])    # SPICE storage
s = s - 35.               # salinity ajusted
T = 1. * np.ones([rw,cl]) # TEMP storage

# (i, j) are (row, colomn) for coefficent matrix
for i in range(6):
    S = np.ones([rw,cl])  # SALT storage
    for j in range(5):
        sp = sp + B(i,j) * T * S
        S  = S * s
    T = T * t

# sp is SPICINESS MATRIX
        
