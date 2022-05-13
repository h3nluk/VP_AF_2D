
import numpy as np
import matplotlib.pyplot as plt

import scipy.signal as sig
import scipy.optimize as opt 


#load data
E_dat = np.genfromtxt('E.dat')
Emax_dat = np.genfromtxt('Emax.dat')
Esq_dat = np.genfromtxt('Esq.dat')

Phi_dat = np.genfromtxt('Phi.dat')
Rho_dat = np.genfromtxt('Rho.dat')

f_dat = np.genfromtxt('f.dat')
f0_dat = np.genfromtxt('f0.dat')

#plotting
fig, axs = plt.subplots(2,3, figsize=(15,15))

#E profile
axs[0,0].plot(E_dat)
axs[0,0].set_title('E(x)')
axs[0,0].grid()

#E_max over time
axs[0,1].plot(Emax_dat[:,0], Emax_dat[:,1])
axs[0,1].set_title(r'$E_{max}(t)$')
axs[0,1].set_yscale('log')
axs[0,1].grid()

#E^2 over over time with Landau Damping fit 
t = Esq_dat[:,0]; E = Esq_dat[:,1]
axs[0,2].plot(t,E)

def func(x,a,b):
	return a*np.exp(b*x)
	
gamma = -0.1533

peaks = sig.find_peaks(E)
popt, pcov = opt.curve_fit(func, t[peaks[0][1:]], E[peaks[0][1:]], p0=[0,gamma])

a = popt[0];  b = popt[1]

axs[0,2].plot(t[peaks[0][1:]], E[peaks[0][1:]], 'r*')
axs[0,2].plot(t, func(t,a,b), color = 'black', linestyle = '--', label= r'$\gamma$ = ' + f'{round(b,4)}')
axs[0,2].grid()
axs[0,2].legend(loc='best')
axs[0,2].set_title(r'$E^2 (t)$')
axs[0,2].set_yscale('log')

#Phi(x)
axs[1,0].plot(Phi_dat)
axs[1,0].set_title(r'$\Phi (x)$')
axs[1,0].grid()

#rho(x)
axs[1,1].plot(Rho_dat)
axs[1,1].set_title(r'$\rho (x)$')
axs[1,1].grid()

#f(x,v) 
# ~ axs[1,2].imshow(np.transpose(f_dat))
# ~ axs[1,2].set_title(r'$f(x,v)$')

#(f_0(x,v)-f(x,v))
axs[1,2].imshow(np.transpose(f0_dat-f_dat),interpolation='gaussian')
axs[1,2].set_title(r'$f_0(x,v)-f(x,v)$')

plt.show()
