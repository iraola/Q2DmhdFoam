'''
USAGE: call the python script adding as an argument the name of the file to read
in this case 'forcesCo/0/forceCoeffs.dat'
COLUMNS TO BE UNPACKED
    runTime
    magUbarStar
    gradPplus
    gradP
    gradPvec (X component)
    repsmaxU
    repsmaxT
'''

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import sys
import os
import pandas as pd

case = os.path.basename(os.getcwd())
fig, ax = plt.subplots(figsize=(12,8))

def animate(i):
    mydata  = pd.read_csv(str(sys.argv[1]), delimiter='\t', header=None, skiprows=1)
    mydata  = mydata.to_numpy()
    time = mydata[:,0]
    Cd = mydata[:,1]
    Cl = mydata[:,2]
    Cm = mydata[:,3]
    ax.clear()
    ax.plot(time, Cl, label='$C_L$')
    # ax.set_title(str(case) +' - $C_L$')
    ax.legend(loc='best')
    ax.set_ylabel('$C_L$')
    ax.set_xlabel('Simulation time (s)')
    ax.grid(linewidth=1, linestyle=':')

ani1 = FuncAnimation(fig, animate, interval=1000)
plt.show()
