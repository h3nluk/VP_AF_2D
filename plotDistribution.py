
import numpy as np
import matplotlib.pyplot as plt

dat = np.genfromtxt('f.dat')
dat0 = np.genfromtxt('f0.dat')

fig, axs = plt.subplots(1,3)

axs[0].imshow(dat)
axs[0].set_title(r'$f(x,v)$')
axs[0].set_xlabel('v')
axs[0].set_ylabel('x')

axs[1].imshow(dat0)
axs[1].set_title(r'$f_0 (x,v)$')
axs[1].set_xlabel('v')
axs[1].set_ylabel('x')

axs[2].imshow(dat-dat0)
axs[2].set_xlabel('v')
axs[2].set_ylabel('x')

plt.show()
