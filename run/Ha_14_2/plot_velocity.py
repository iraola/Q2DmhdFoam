import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.animation as animation
import sys
import os
import pandas as pd

font = {'size' : 10,
        'family': 'serif',
        }
matplotlib.rc('font', **font)


case= os.path.basename(os.getcwd())
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

D   = 0.2
v0  = 0.00021866

#def animate(i):
ct, vy1, vy2, vy3  = np.loadtxt('postProcessing/probes/450000/U', unpack=True, skiprows=5, usecols = (0,3,6,9), dtype=str)
ct  = ct.astype(np.float)
vy1 = vy1.astype(np.float)
vy2 = vy2.astype(np.float)
vy3 = vy3.astype(np.float)
#=============================================================================================
#ax1.clear()

ax1.plot(  ct/(D/v0), vy1/v0, '-', label='0.2D')
ax1.plot(  ct/(D/v0), vy2/v0, '-', label='0.7D', color='black')
ax1.plot(  ct/(D/v0), vy3/v0, '-', label='1.2D', color='green')

ax1.set_title(' Vertical velocity suppression - Transition to Q2D Flow')
ax1.legend(loc='lower right')
ax1.set_xlabel('dimensionless time '+r'$(t/(D/U)) $') #('time (s)')
ax1.set_ylabel('dimensionless vertical velocity '+r'$(v_y/U)$')
ax1.grid(linewidth=1, linestyle=':')
ax1.set_ylim(top=0.26,bottom=-0.15)
#ax1.set_xlim(left=0,right=95/(D/v0))
#ax1.set_xlim(left=400,right=95/(D/v0))
ax1.set_xlim(left=400, right=560)
# -------------------------------------------------------------------

fig.tight_layout()
plt.savefig('mhdtransition.png', format='png', dpi=1200)
#ani1 = animation.FuncAnimation(fig, animate, interval=1000)
plt.show()
