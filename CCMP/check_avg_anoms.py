import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as pltd

def running_mean(x,N):
    
    cumsum = np.cumsum(np.insert(x, 0, 0)) 

    return (cumsum[N:] - cumsum[:-N]) / N 

###################################
# f = open('vwind_anomaly_stats','r')
f = open('vwind_stats','r')
A = f.read()
A2 = A.split() 
A3 = [float(x) for x in A2]
A3 = A3[5:]
B = np.reshape(A3,(-1,6))

# f2 = open('uwind_anomaly_stats','r')
f2 = open('uwind_stats','r')
A = f2.read()
A2 = A.split()
A3 = [float(x) for x in A2]
A3 = A3[5:]
C = np.reshape(A3,(-1,6))

# YEAR, MON, DAY, MEAN, STD, MIN, MAX

wind = running_mean(B[:,3],365)
wind2 = running_mean(C[:,3],365)

#print '365d vwind minimum begins:', B[np.where(wind==np.min(wind))[0],:3]
#print '365d uwind minimum begins:', C[np.where(wind2==np.min(wind2))[0],:3]

#print np.where(wind==np.min(wind))
#print np.where(wind2==np.min(wind2))

#print wind[1450:1470]
#print wind2[1450:1470]

date_vals = np.zeros(B.shape[0])
for nt in range(B.shape[0]):
    day_time = dt.datetime(int(B[nt,0]),1,1) + dt.timedelta(days = int(B[nt,1]))
    date_vals[nt] = pltd.date2num(day_time)


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(date_vals[:len(wind)],wind)
ax.plot([date_vals[0],date_vals[len(wind)-1]],np.mean(wind)*np.ones(2))
ax.plot(date_vals[:len(wind2)],wind2)
ax.plot([date_vals[0],date_vals[len(wind2)-1]],np.mean(wind2)*np.ones(2))
ax.xaxis_date()


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx() 

plt.title('Daily averaged zonal (cyan) and meridional (blue) CCMP anomalies + ERAi Climatologies winds (m/s) over CCS domain \n Black lines indicate monthly means')

fig2 = plt.figure()
ax3 = fig2.add_subplot(111)
ax4 = ax3.twinx()

ax1.plot(date_vals,B[:,2]) #vwind
ax2.plot(date_vals,C[:,2],'-c') #uwind

ax1.set_ylim(-22,8)
ax2.set_ylim(-6,22)

ax3.set_ylim(-22,8)
ax4.set_ylim(-6,22)
nmdays = [31,28,31,30,31,30,31,31,30,31,30,31]
d0 = 0
yd0 = 0
#print len(date_vals)
for ny in range(1990,2010+1):
    #if ny % 4 == 0:
    #   nmdays[1]=29
    #else:
    #   nmdays[1]=28
    ydn = yd0 + 365 
    for nm in range(12):
        dn = d0 + nmdays[nm]
        ax1.plot([date_vals[d0],date_vals[dn-1]],np.mean(B[d0:dn,2])*np.ones(2),'-k')
        ax2.plot([date_vals[d0],date_vals[dn-1]],np.mean(C[d0:dn,2])*np.ones(2),'-k')
        d0 = dn

    ax3.plot([date_vals[yd0],date_vals[ydn-1]],np.mean(B[yd0:ydn,2])*np.ones(2),'-k')
    ax4.plot([date_vals[yd0],date_vals[ydn-1]],np.mean(C[yd0:ydn,2])*np.ones(2),'-k') 

    Berr = np.std(B[yd0:ydn,2])
    Cerr = np.std(C[yd0:ydn,2])

    ax3.fill_between([date_vals[yd0],date_vals[ydn-1]], np.mean(B[yd0:ydn,2])-Berr, np.mean(B[yd0:ydn,2])+Berr,facecolor='lightskyblue')
    ax4.fill_between([date_vals[yd0],date_vals[ydn-1]], np.mean(C[yd0:ydn,2])-Cerr, np.mean(C[yd0:ydn,2])+Cerr,facecolor='darkorange')

    yd0 = ydn    

ax1.plot([date_vals[0],date_vals[-1]],np.mean(B[:,2])*np.ones(2),'--g')
ax2.plot([date_vals[0],date_vals[-1]],np.mean(C[:,2])*np.ones(2),'--g')
ax3.plot([date_vals[0],date_vals[-1]],np.mean(B[:,2])*np.ones(2),'--r')
ax4.plot([date_vals[0],date_vals[-1]],np.mean(C[:,2])*np.ones(2),'--r')
#ax1.set_yticks([2,4,6])
#ax2.set_yticks([2,4,6])
ax1.xaxis_date()
ax3.xaxis_date()
plt.show()
