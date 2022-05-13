
import numpy as np
import matplotlib.pyplot as plt

dat = np.genfromtxt('f.dat')
dat0 = np.genfromtxt('f0.dat')
datdiff = dat-dat0

sizex, sizev = np.shape(dat)

fig, axs = plt.subplots(2,3)

axs[0,0].imshow(dat)
axs[0,0].axhline(sizex//2, color='blue')
axs[0,0].axvline(sizev//2, color='red')
axs[0,0].set_title(r'$f(x,v)$')

axs[0,1].imshow(dat0)
axs[0,1].axhline(sizex//2, color='blue')
axs[0,1].axvline(sizev//2, color='red')
axs[0,1].set_title(r'$f_0(x,v)$')

axs[0,2].imshow(datdiff)
axs[0,2].axhline(sizex//2, color='blue')
axs[0,2].axvline(sizev//2, color='red')
axs[0,2].set_title(r'$f - f_0$')

axs[1,0].plot(dat[:, sizev//2], color = 'red')
axs[1,0].plot(dat[sizex//2, :], color = 'blue')

axs[1,1].plot(dat0[:, sizev//2], color = 'red')
axs[1,1].plot(dat0[sizex//2, :], color = 'blue')

axs[1,2].plot(datdiff[:, sizev//2], color = 'red')
axs[1,2].plot(datdiff[sizex//2, :], color = 'blue')

plt.show()
