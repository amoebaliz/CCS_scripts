import numpy as np
import matplotlib.pyplot as plt



A = np.load('/Users/elizabethdrenkard/Desktop/LENS_clim_difs.npy')

for nl in range(A.shape[1]):
    if nl == (34-1):
       colr = 'red'; zord = 4
    elif nl == (19-1):
       colr = 'blue'; zord = 4
    elif nl == (11-1):
       colr = 'green'; zord = 4 
    elif nl == (17-1):
       colr = 'black'; zord = 4
    elif nl == (18-1):
       colr = 'purple'; zord = 4

    elif nl == (24-1):
       colr = 'pink'; zord = 4
    else:
       colr = 'grey'; zord = 2


    plt.plot(np.arange(1,12+1),A[:,nl],color=colr,marker='o',zorder=zord)


plt.figure()
plt.plot(np.arange(1,12+1),A[:,15-1],color='blue',marker='o',zorder=2)
plt.plot(np.arange(1,12+1),A[:,34-1],color='green',marker='o',zorder=2)
plt.plot(np.arange(1,12+1),A[:,24-1],color='purple',marker='o',zorder=2)
plt.plot(np.arange(1,12+1),A[:,17-1],color='olive',marker='o',zorder=2)
plt.plot(np.arange(1,12+1),A[:,13-1],color='red',marker='o',zorder=2)
plt.show()

