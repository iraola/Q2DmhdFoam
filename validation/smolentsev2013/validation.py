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

# Config plots
font_config = {'font.family'     : 'serif',
               'font.size'       : 16,
               'lines.linewidth' : 2}
plt.rcParams.update(font_config)

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
    'a'   : 2/2,
    'b'   : 2/2,
    'q0'  : 0,
    'qWall' : 0,
    # CAREFUL WITH 'm': m(Q2D) = dimless. m(tepot) = meters
    #                   therefore m(Q2D) = m(tepot) * characteristic length
    'm'   : 6.3 * 0.15/2,
    'Th'  : '?',
    'Tc'  : '?',
    'g'   : -9.81}
phys_dict = {
    'rho0'   : 9720,
    'nu'     : 1.54012345679012e-07,
    'Cp'     : 189,
    'k'      : 15.14,
    'beta'   : 1.2e-4,
    'sigma'  : 763000}
mesh_dict = {       # SETUP FOR CYCLIC ONE-CELL CASE
    'Lx'    : 60,
    'LxHalf': 60/2,
    'Nx'    : 100,
    'c2c_bl'    : 1.05,
    'c2c_bulk'  : 1.001
}
a = tag_dict['a']   # take 'a' length that we'll use later

# Initialize plot
fig_scaled, ax_scaled = plt.subplots(figsize=(12,6))
fig, ax = plt.subplots(figsize=(12,6))

"""
def getConditions(filename):
    ''' Get Ha, Re, Gr from the name of a file as real numbers'''
    my_list = filename.split('_')
    Ha = float(my_list[3])
    Re = float(my_list[4])
    Gr = float(my_list[5])
    return Ha, Re, Gr
"""
def getConditions():
    return 50, 5000, 1e8

### LOOP
# Remove old directories present
sp.call("rm -r -f case_*", shell=True)
errors = {}
i = 0
for file in os.listdir(validation_dir):
    # Get dimless numbers from the filenames
    Ha, Re, Gr = getConditions()
    print('\nHa={}, Re={}, Gr={}'.format(Ha, Re, Gr))

    # Run case with the meshAndGo module
    meshAndGo(Ha, Re, Gr, mesh_dict=mesh_dict, tag_dict=tag_dict,
	      phys_dict=phys_dict)

    # Get latest time directory name and load simulation data
    latest_time = getLatestTime(postprocess_dir)
    filename = postprocess_dir + latest_time + '/' + postprocess_file
    z, U, _, _ = np.loadtxt(filename, unpack= True)
    z_val, _, U_val = np.loadtxt(validation_dir + '/' + file, unpack= True)
    label_q2d = f'Ha = {Ha:.0f} - Q2D'
    label_val = f'Ha = {Ha:.0f} - Analytical'

    # Plot unscaled data
    ax.plot(z, U, linestyle='-', color=color_q2d[i], label=label_q2d)
    ax.plot(z_val, U_val, linestyle='--', color=color_val[i], label=label_val)

    # Normalize data
    #  for simulation data
    U_scaled = U / U.mean()
    z_scaled = (z - a) / a
    #  for validation data
    U_val_scaled = U_val / U_val.mean()
    z_val_scaled = (z_val - a) / a
    # Plot scaled data
    ax_scaled.plot(z_scaled, U_scaled, linestyle='-', color=color_q2d[i], label=label_q2d)
    ax_scaled.plot(z_val_scaled, U_val_scaled, linestyle='--', color=color_val[i], label=label_val)

    # Get metric of performance:
    #   root mean squared error over the normalized velocities
    rmse = np.sqrt(mean_squared_error(y_true=U_val_scaled, y_pred=U_scaled))
    print('Root mean squared error of the normalized velocity is', rmse)
    errors[i+1] = [Ha, Re, Gr, rmse]
    # Store the simulated case in other directory
    sp.call("mv case case_" + str(Ha) + "_" + str(Gr), shell=True)
    # Update counter
    i += 1

# Write validation errors file
pd.DataFrame.from_dict(errors, orient='index',
                       columns=['Ha', 'Re', 'Gr', 'rmse']).to_csv(
                       'validation_error_tepot.csv', index=False)

# Plot configuration
ax_scaled.set_ylabel('Normalized velocity, $U/\overline{U}$')
ax_scaled.set_xlabel('Dimensionless channel length, $y/a$')
ax_scaled.grid(True)
ax_scaled.legend(loc='best')
fig_scaled.savefig('validation_tepot_scaled.png', format='png',dpi=300)

ax.set_ylabel('Velocity (m/s)')
ax.set_xlabel('Channel length (m)')
ax.grid(True)
ax.legend(loc='best')
fig.savefig('validation_tepot.png', format='png',dpi=300)
plt.show()