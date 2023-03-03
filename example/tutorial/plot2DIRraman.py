import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

Data = np.loadtxt('2DIRraman.I.par.dat')

w1Dim=int(np.sqrt(len(Data[:,0])))
w3Dim=w1Dim

wmin=1125
wmax=1275

step=25
mmax=np.max(np.abs(Data[:,3]))
Data[:,3]=Data[:,3]/mmax
Data[1,3]=1
Data[2,3]=-1

w1 = np.reshape(Data[:,0],(w1Dim,w3Dim))
w3 = np.reshape(Data[:,1],(w1Dim,w3Dim))

PlotMap = np.reshape(Data[:,3],(w1Dim,w3Dim))

plt.contourf(w1,w3,PlotMap,20,cmap=plt.cm.bwr)
plt.xlabel('$\omega_1$ [cm$^{-1}$]',fontsize=16)
plt.ylabel('$\omega_3$ [cm$^{-1}$]',fontsize=16)
plt.gca().set_aspect('equal')

plt.xlim(wmin,wmax)
plt.ylim(wmin,wmax)

plt.plot([ wmin , wmax],[wmin,wmax],color='black')
plt.xticks(np.arange(wmin,wmax+1,step),fontsize=12,rotation=0)
plt.yticks(np.arange(wmin,wmax+1,step),fontsize=12,rotation=0)
plt.show()

Data = np.loadtxt('2DIRraman.II.par.dat')

w1Dim=int(np.sqrt(len(Data[:,0])))
w3Dim=w1Dim

wmin=1125
wmax=1275
wmin1=wmin+1200
wmax1=wmax+1200

step=25
mmax=np.max(np.abs(Data[:,3]))
Data[:,3]=Data[:,3]/mmax
Data[1,3]=1
Data[2,3]=-1

w1 = np.reshape(Data[:,0],(w1Dim,w3Dim))
w3 = np.reshape(Data[:,1],(w1Dim,w3Dim))

PlotMap = np.reshape(Data[:,3],(w1Dim,w3Dim))

plt.contourf(w1,w3,PlotMap,20,cmap=plt.cm.bwr)
plt.xlabel('$\omega_1$ [cm$^{-1}$]',fontsize=16)
plt.ylabel('$\omega_3$ [cm$^{-1}$]',fontsize=16)
plt.gca().set_aspect('equal')

plt.xlim(wmin1,wmax1)
plt.ylim(wmin,wmax)

plt.plot([ wmin1 , wmax1],[wmin,wmax],color='black')
plt.xticks(np.arange(wmin1,wmax1+1,step),fontsize=12,rotation=0)
plt.yticks(np.arange(wmin,wmax+1,step),fontsize=12,rotation=0)
plt.show()

