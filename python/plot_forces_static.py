'''
USAGE: call the python script adding as an argument the name of the file to read
in this case 'forcesCo/0/forceCoeffs.dat'
'''

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import sys
import os
import pandas as pd

case_dir = os.getcwd()
fig, ax = plt.subplots(figsize=(12,8))
forces_dir = os.path.join(case_dir, sys.argv[1])
data = pd.DataFrame()
for subdir in os.listdir(forces_dir):
    print(subdir, forces_dir)
    time_dir = os.path.join(forces_dir, subdir)
    print(time_dir)
    filename = os.listdir(time_dir)[0]
    file_path = os.path.join(time_dir, filename)
    # Read file
    data = pd.concat([data, pd.read_csv(file_path, delimiter='\t', header=None, skiprows=1)],
                     axis=0)

data  = data.to_numpy()
time = data[:,0]
Cd = data[:,1]
Cl = data[:,2]
Cm = data[:,3]
ax.clear()
ax.plot(time, Cl, label='$C_L$')
# ax.set_title(str(case) +' - $C_L$')
ax.legend(loc='best')
ax.set_ylabel('$C_L$')
ax.set_xlabel('Simulation time (s)')
ax.grid(linewidth=1, linestyle=':')

plt.show()
