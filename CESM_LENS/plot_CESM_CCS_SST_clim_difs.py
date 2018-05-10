import numpy as np
import matplotlib.pyplot as plt



A = np.load('/Users/elizabethdrenkard/Desktop/LENS_clim_difs.npy')

avg_A = np.mean(A,axis=0)
print np.argsort(avg_A)+1
print avg_A[np.argsort(avg_A)]

for nl in range(A.shape[1]):
    if nl == (17-1):
       colr = 'red'; mec_clr = 'firebrick'; zord = 4
    elif nl == (19-1):
       colr = 'blue'; mec_clr = 'darkblue'; zord = 4

    # COMPARITIVE COLOR
    # elif nl == (34-1):
    #   colr = 'cyan'; zord = 4 

    else:
       colr = 'silver'; mec_clr = 'darkgrey'; zord = 2


    plt.plot(np.arange(1,12+1),A[:,nl],marker='o', ms=10, mfc=colr,mec=mec_clr,ls='None',zorder=zord)


labels = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG','SEP','OCT','NOV','DEC']
plt.xticks(np.arange(1,12+1),labels)
plt.show()

