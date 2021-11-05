import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import sys
import os
import pandas as pd

case= os.path.basename(os.getcwd())
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

D   = 0.2
v0  = 1
ax2 = ax1.twinx()
def animate(i):
    ct, cd, cl = np.loadtxt('forces/0/forceCoeffs.dat', delimiter='\t', unpack=True, skiprows=5, usecols = (0,1,2))
#    deltaT   = [t-s for s, t in zip(x, x[1:])]
# utilitza PANDAS (teoricament mes rapid que numpy (np) =========================================
#    mydata    = pd.read_csv('forces/0/forceCoeffs.dat', header=None, delimiter='\t',  skiprows=5, usecols = [0,1,2])
#    ct        =   mydata[0]
#    cd        =   mydata[1]
#    cl        =   mydata[2]
#    print(ct[:10])
#    print(cd[:10])
#=============================================================================================

    ax1.clear()
    ax2.clear()
#    ax3.clear()
#    ax4.clear()

    ax2.plot(  ct/(D/v0), 4*cd, '-', label='Cd')  # aplico dimensionless time i el cd l'aplico només a la cara d'entrada del cilindre, llavors x4
    ax1.plot(  ct/(D/v0), 4*cl, '-', label='Cl', color='black')  # aplico dimensionless time i el cd l'aplico només a la cara d'entrada del cilindre, llavors x4

    ax1.set_title(' Drag and Lift Coefficients - Hydrodynamic Flow')
#    ax3.legend()
    ax1.legend(loc='center right')
    ax2.legend(loc='lower right')
    ax2.set_ylabel('Cd')
    ax1.set_ylabel('Cl')
    ax1.set_xlabel('dimensionless time (t/(D/v0) ') #('time (s)')
    ax1.grid(linewidth=1, linestyle=':')
#    ax2.set_yscale('log')
    ax2.set_ylim(top=2,bottom=0)
    ax1.set_ylim(top=1.25,bottom=-0.75)
# -------------------------------------------------------------------

ani1 = animation.FuncAnimation(fig, animate, interval=1000)
plt.show()
