#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import subprocess as sp
import pdb
# import custom modules
sys.path.insert(1, '../../python')
from meshAndGo import meshAndGo
from MHDutils import getLatestTime

### INITIALIZATIONS
validation_file = 'vetcha2009_Ha50_Re1e4.csv'
postprocess_dir = 'case/postProcessing/sets/'
postprocess_file = 'line_centreProfile_U.xy'
# Case setup tags
tag_dict = {
    'B'   : '?',
    'Ux'  : '?',
    'T0'  : 0.0,
    'a'   : 0.15/2,
    'b'   : 0.15/2,
    'q0'  : 0,
    'qWall' : 0,
    'm'   : 1,
    'Th'  : '?',
    'Tc'  : '?',
    'g'   : -9.81}
phys_dict = {
    'rho0'   : 9720,
    'nu'     : 1.54e-7,
    'Cp'     : 189,
    'k'      : 22.36,
    'beta'   : 1.2e-4,
    'sigma'  : 763000}
mesh_dict = {       # SETUP FOR CYCLIC ONE-CELL CASE
    'Lx'    : 0.1,
    'LxHalf': 0.1/2,
    'Nx'    : 1}
Grashofs = [5e6, 5e7, 1e8]
Re = 1e4
Ha = 50
# Load validation dataset
names = ['x_Gr5e6','y_Gr5e6','x_Gr5e7','y_Gr5e7','x_Gr1e8','y_Gr1e8']
val_data = pd.read_csv(validation_file, skiprows=2, names=names)
val_data.head()
# Initialize plot
fig, ax = plt.subplots(figsize=(12,6))

### LOOP
# Remove old directories present
sp.call("rm -r -f case_*", shell=True)
for Gr in Grashofs:
    i = Grashofs.index(Gr)
    # Run case with the meshAndGo module
    meshAndGo(Ha, Re, Gr,
        mesh_dict=mesh_dict, tag_dict=tag_dict, phys_dict=phys_dict)
    # Store the simulated case in other directory
    sp.call("cp -r case case_" + str(Gr), shell=True)
    # Get latest time directory name and load simulation data
    latest_time = getLatestTime(postprocess_dir)
    filename = postprocess_dir + latest_time + '/' + postprocess_file
    # Plot sim. data
    z, U, _, _ = np.loadtxt(filename, unpack= True)
    # Normalize data
    U /= U.mean()
    z = (z - tag_dict['a']) / tag_dict['a']
    # Plot simulationdata
    ax.plot(z, U, label='Q2D Gr='+str(Gr))
    # Plot validation data
    ax.plot(val_data[val_data.columns[2*i]],
            val_data[val_data.columns[2*i+1]], label='Vetcha Gr='+str(Gr))


ax.legend(loc='best')
plt.show()
