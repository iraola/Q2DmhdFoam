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

# tepot check
tepot_files = ['lineSampled_theta_Ux_250_500_0',
               'lineSampled_theta_Ux_1000_500_0',
               'lineSampled_theta_Ux_2500_500_0']

# FUNCTIONS
# def shercliff_profile(Ha, gradP, nu, a, b, z):
#     ''' Mas de les Valls citation of Ni et al. (2007) corrected solution '''
#
#     print('\nRunning analytical procedure to calculate Shercliff profile',
#            'for the validation case:\n',
#            'Ha = {}, gradP = {}, nu = {}, a = {}, b = {}, z = {} to {}\n'
#            .format(Ha, gradP, nu, a, b, min(z), max(z)))
#
#     y = 0       # center of the channel, study side boundary layer
#     # Initialization
#     eta = y / a     # dimensionless coordinate along magnetic field lines
#     chi = z / a     # dimensionless coordinate perpendicular to magnetic field lines
#     l = b / a       # aspect ratio
#     db = 0.0        # Hartmann wall conductivity ratio, zero for Shercliff case
#     tol = 1e-21
#     # Loop
#     V = np.zeros(z.shape)
#     V_old = V + 1
#     k = 0
#     iter = 0
#     while abs((V_old - V).mean() / V_old.mean()) > tol:
#         iter += 1
#         V_old = V.copy()
#         alpha = (k + 1 / 2) * np.pi / l
#         r1 = 1 / 2 * (Ha + (Ha**2 + 4 * alpha**2) ** (1 / 2))
#         r2 = 1 / 2 * (-Ha + (Ha**2 + 4 * alpha**2) ** (1 / 2))
#         N = (Ha**2 + 4 * alpha**2) ** (1 / 2)
#         V2 = ((db * r2 + (1 - exp(-2 * r2)) / (1 + exp(-2 * r2)))
#               * (exp(-r1 * (1 - eta)) + exp(-r1 * (1 + eta))) / 2) \
#              / ((1 + exp(-2 * r1)) / 2 * db * N
#                 + (1 + exp(-2 * (r1 + r2))) / (1 + exp(-2 * r2)))
#         V3 = ((db * r1 + (1 - exp(-2 * r1)) / (1 + exp(-2 * r1)))
#               * (exp(-r2 * (1 - eta)) + exp(-r2 * (1 + eta))) / 2) \
#              / ((1 + exp(-2 * r2)) / 2 * db * N
#                 + (1 + exp(-2 * (r1 + r2))) / (1 + exp(-2 * r1)))
#         V += (2 * (-1)**k * np.cos(alpha * chi)) / (l * alpha**3) * (1 - V2 - V3)
#         k += 1
#     print('Ran {} iterations'.format(iter))
#     U = nu**-1 * V * (-gradP) * a**2
#     if U.mean() < 0: U = -U
#     return U

# ####################### test function
# U_0 = 0.0011
# Ha = 5000
# gradP = 1e-5
# z_val = np.linspace(-0.075, 0.075, 1000) # 'y' in Elisabet's
# U_val = shercliff_profile(Ha, gradP, nu, a, b, z_val)
# ax_scaled.plot(z_val, U_val, linestyle='--')
# print(max(abs(U_val)))
# #######################

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

    # # Load last value of pressure gradient
    # _, _, _, gradP, _, _ = np.loadtxt('case/output.out', unpack=True)
    # gradP = gradP[-1]
    # # Calculate validation data with Shercliff procedure
    # z_val = z - a   # center axis data around zero to run shercliff
    # U_val = shercliff_profile(Ha, gradP, nu, a, b, z_val)
    # # Plot data without normalization
    # plt.figure()
    # plt.plot(z, U_val, label='validation')
    # plt.plot(z, U, label='simulation')
    # plt.legend(loc='best')
    # plt.title('Hartmann ' + str(Ha))

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
