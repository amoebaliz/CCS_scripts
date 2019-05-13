import numpy as np
import matplotlib.pyplot as plt



A = np.ones([4,4])
A = A[1,:]*2
A = A[2,:]*3
A = A[3,:]*4


lon = np.tile([1,2,3,4],(4,1))
lat = np.tile([[1],[2],[3],[4]],(1,4))


