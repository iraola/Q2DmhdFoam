import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import sys
import os
import pandas as pd

case = os.path.basename(os.getcwd())
fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)

'''
COLUMNS TO BE UNPACKED
    runTime
    magUbarStar
    gradPplus
    gradP
    gradPvec (X component)
    repsmaxU
    repsmaxT
'''

def animate(i):
    mydata  = pd.read_csv(str(sys.argv[1])+'/output.out', delimiter='\t',
        header=None, skiprows=5)
    mydata  = mydata.to_numpy()
    runTime         =   mydata[:,0]
    magUbarStar     =   mydata[:,1]
    gradPplus       =   mydata[:,2]
    gradP           =   mydata[:,3]
    repsmaxU        =   mydata[:,4]
    repsmaxT        =   mydata[:,5]
#==============================================================================
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()

    ax1.plot(runTime, magUbarStar, '-', label='magUbarStar')
    ax2.plot(runTime, gradPplus, '-', label='gradPplus'  , color='r')
    ax3.plot(runTime, gradP, '-', label='gradP', color='k')
    ax4.plot(runTime, repsmaxU, '-', label='repsmaxU', color='b')
    ax4.plot(runTime, repsmaxT, '-', label='repsmaxT', color='r')

    #ax1.set_title(str(case) +' - output.out')
    fig.suptitle(str(case) +' - output.out', fontsize=16)
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    ax3.legend(loc='best')
    ax4.legend(loc='best')
    ax1.set_ylabel('magUbar')
    ax2.set_ylabel('gradP')
    ax3.set_ylabel('gradP')
    ax4.set_ylabel('reps max')
    ax1.set_xlabel('Simulation time (s)') #('time (s)')
    ax2.set_xlabel('Simulation time (s)') #('time (s)')
    ax3.set_xlabel('Simulation time (s)') #('time (s)')
    ax4.set_xlabel('Simulation time (s)') #('time (s)')

    ax1.grid(linewidth=1, linestyle=':')
    ax2.grid(linewidth=1, linestyle=':')
    ax3.grid(linewidth=1, linestyle=':')
    ax4.grid(linewidth=1, linestyle=':')
#    ax3.grid(linewidth=1, linestyle=':')
    ax4.set_yscale('log')
#    ax1.set_ylim(top=10,bottom=-1)
#    ax3.set_ylim(top=0.04,bottom=0)
# -----------------------------------------------------------------------------

ani1 = FuncAnimation(fig, animate, interval=1000)
plt.show()
