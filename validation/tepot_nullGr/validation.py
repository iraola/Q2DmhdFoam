#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import subprocess as sp
from sklearn.metrics import mean_squared_error
import pdb
# import custom modules
sys.path.insert(1, '../../python')
from meshAndGo import meshAndGo
from MHDutils import getLatestTime

### INITIALIZATIONS
validation_dir = 'samples'
postprocess_dir = 'case/postProcessing/sets/'
postprocess_file = 'line_centreProfile_U.xy'
color_q2d = ['b', 'g', 'r', 'm', 'darkorange']                        # strong colors
color_val = ['royalblue', 'limegreen', 'salmon', 'violet', 'orange']  # light colors
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
    'nu'     : 1.5401234567901237e-07,
    'Cp'     : 189,
    'k'      : 15.14,
    'beta'   : 1.2e-4,
    'sigma'  : 763000}
mesh_dict = {       # SETUP FOR CYCLIC ONE-CELL CASE
    'Lx'        : 0.1,
    'LxHalf'    : 0.1/2,
    'Nx'        : 1,
    'c2c_bl'    : 1.05,
    'c2c_bulk'  : 1.001
}
a = tag_dict['a']   # take 'a' length that we'll use later

# Initialize plot
fig, ax = plt.subplots(figsize=(12,6))

def getConditions(filename):
    ''' Get Ha, Re, Gr from the name of a file as real numbers'''
    my_list = filename.split('_')
    Ha = float(my_list[3])
    Re = float(my_list[4])
    Gr = float(my_list[5])
    return Ha, Re, Gr

### LOOP
# Remove old directories present
sp.call("rm -r -f case_*", shell=True)
errors = {}
i = 0
for file in os.listdir(validation_dir):
    # Get dimless numbers from the filenames
    Ha, Re, Gr = getConditions(file)
    print('Ha={}, Re={}, Gr={}'.format(Ha, Re, Gr))
    # Run case with the meshAndGo module
    meshAndGo(Ha, Re, Gr,
        mesh_dict=mesh_dict, tag_dict=tag_dict, phys_dict=phys_dict)
    # Get latest time directory name and load simulation data
    latest_time = getLatestTime(postprocess_dir)
    filename = postprocess_dir + latest_time + '/' + postprocess_file
    # Plot
    z, U, _, _ = np.loadtxt(filename, unpack= True)
    # Conditional to allow input sample files with different column formats
    # (mainly with or without 'theta' field)
    if np.loadtxt(validation_dir + '/' + file).shape[1] > 2:
        z_val, _, U_val = np.loadtxt(validation_dir + '/' + file, unpack= True)
    else:
        z_val, U_val = np.loadtxt(validation_dir + '/' + file, unpack= True)
    # Normalize data
    # for simulation data
    U /= U.mean()
    z = (z - a) / a
    # for validation data
    U_val /= U_val.mean()
    z_val = (z_val - a) / a
    # Plot simulation data
    label_q2d = 'Q2D Ha=' + str(Ha) + ' Gr='+str(Gr)
    ax.plot(z, U, linestyle='-', color=color_q2d[i], label=label_q2d)
    # Plot validation data
    label_val = 'tepot Ha=' + str(Ha) + ' Gr='+str(Gr)
    ax.plot(z_val, U_val, linestyle='--', color=color_val[i], label=label_val)
    # Get metric of performance:
    #   root mean squared error over the normalized velocities
    rmse = np.sqrt(mean_squared_error(y_true=U_val, y_pred=U))
    print('Root mean squared error of the normalized velocity is', rmse)
    errors[i+1] = [Ha, Re, Gr, rmse]
    # Store the simulated case in other directory
    sp.call("mv case case_" + str(Ha) + "_" + str(Gr), shell=True)
    # Update counter
    i += 1

# Write validation errors file
pd.DataFrame.from_dict(errors, orient='index',
                       columns=['Ha', 'Re', 'Gr', 'rmse']).to_csv(
                       'validation_error_shercliff.csv', index=False)

ax.legend(loc='best')
fig.savefig('validation_tepot_nullGr.png',format='png',dpi=200)
plt.show()
