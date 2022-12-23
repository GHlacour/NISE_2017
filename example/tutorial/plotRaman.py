import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

Data = np.loadtxt('TD_Raman_VV.dat')
Data2 = np.loadtxt('TD_Raman_VH.dat')

plt.plot(Data[:,0],Data[:,1],label='VV')
plt.plot(Data[:,0],Data2[:,1],label='VH')
plt.xlabel('Time [fs]',fontsize=16)
plt.ylabel('Response function [arb.u.]',fontsize=16)
plt.legend()
plt.show()

Data = np.loadtxt('Raman_VV.dat')
Data2 = np.loadtxt('Raman_VH.dat')

plt.plot(Data[:,0],Data[:,1],label='VV')
plt.plot(Data[:,0],Data2[:,1],label='VH')
plt.xlabel('$\omega$ [cm$^{-1}$]',fontsize=16)
plt.ylabel('Absorption [arb.u.]',fontsize=16)
plt.legend()
plt.show()
