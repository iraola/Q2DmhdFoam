#!/usr/bin/env python3
import math
import numpy as np
from numpy import exp
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import subprocess as sp
from sklearn.metrics import mean_squared_error
from pdb import set_trace
# import custom modules
sys.path.insert(1, '../../python')
from meshAndGo import meshAndGo
from MHDutils import getLatestTime, shercliff_profile

### INITIALIZATIONS
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
# Do not provide mesh_dict. blockMeshDict has no tags
a = tag_dict['a']
b = tag_dict['b']
nu = phys_dict['nu']

# Initialize plot
fig_scaled, ax_scaled = plt.subplots(figsize=(12,6))
fig, ax = plt.subplots(figsize=(12,6))

# Initialize conditions: Q2D laminar conditions according to Smolentsev
Ha_list = [250, 1000, 2500]
Re = 500
Gr = 0.0

### LOOP
# Remove old directories present
sp.call("rm -r -f case_*", shell=True)
errors = {}
i = 0
for Ha in Ha_list:
    print('\nHa={}, Re={}, Gr={}'.format(Ha, Re, Gr))
    # Run case with the meshAndGo module
    meshAndGo(Ha, Re, Gr, mesh_dict=mesh_dict, tag_dict=tag_dict,
              phys_dict=phys_dict)
    # Get latest time directory name and load simulation data
    latest_time = getLatestTime(postprocess_dir)
    filename = postprocess_dir + latest_time + '/' + postprocess_file
    # Load simulation data
    z, U, _, _ = np.loadtxt(filename, unpack=True)

    # Plot validation versus simulation data one per run
    z_val, U_val = np.loadtxt('samples/' + tepot_files[i], unpack=True)
    plt.figure()
    plt.plot(z, U, label='Q2D')
    plt.plot(z_val, U_val, label='validation')
    plt.legend(loc='best')
    plt.title('Hartmann ' + str(Ha))

    # Show absolute values before normalization
    print('Simulation max. velocity is', U.max())
    print('Validation max. velocity is', U_val.max())

    # General big plotting
    # Prepare labels
    label_q2d = 'Q2D Ha=' + str(Ha) + ' Gr='+str(Gr)
    label_val = 'Analytical Ha=' + str(Ha) + ' Gr='+str(Gr)
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
                       'validation_error_shercliff.csv', index=False)

# Plot configuration
ax_scaled.set_ylabel('Normalized velocity, $U/\overline{U}$')
ax_scaled.set_xlabel('Dimensionless channel length, $y/a$')
ax_scaled.grid(True)
ax_scaled.legend(loc='best')
fig_scaled.savefig('validation_shercliff_scaled.png', format='png',dpi=300)

ax.set_ylabel('Velocity (m/s)')
ax.set_xlabel('Channel length (m)')
ax.grid(True)
ax.legend(loc='best')
fig.savefig('validation_shercliff.png', format='png',dpi=300)
plt.show()
sp.call("python3 IDM_Postprocess.py")
