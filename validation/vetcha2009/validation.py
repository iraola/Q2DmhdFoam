#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import subprocess as sp
from sklearn.metrics import mean_squared_error
import pdb
# import custom modules
sys.path.insert(1, '../../python')
from meshAndGo import meshAndGo
from MHDutils import getLatestTime

### INITIALIZATIONS
validation_file = 'vetcha2009_Ha50_Re1e4.csv'
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
    'nu'     : 1.54012345679012e-07,
    'Cp'     : 189,
    'k'      : 22.36,
    'beta'   : 1.2e-4,
    'sigma'  : 763000}
mesh_dict = {       # SETUP FOR CYCLIC ONE-CELL CASE
    'Lx'    : 0.1,
    'LxHalf': 0.1/2,
    'Nx'    : 1,
    'c2c_bl'    : 1.05,
    'c2c_bulk'  : 1.001
}
a = tag_dict['a']   # take 'a' length that we'll use later
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
errors = {}
i = 0
for Gr in Grashofs:
    i = Grashofs.index(Gr)
    # Run case with the meshAndGo module
    meshAndGo(Ha, Re, Gr,
        mesh_dict=mesh_dict, tag_dict=tag_dict, phys_dict=phys_dict)
    # Get latest time directory name and load simulation data
    latest_time = getLatestTime(postprocess_dir)
    filename = postprocess_dir + latest_time + '/' + postprocess_file
    # Plot
    z, U, _, _ = np.loadtxt(filename, unpack= True)
    # Normalize data
    U /= U.mean()
    z = (z - a) / a
    # Plot simulation data
    label_q2d = 'Q2D '+ 'Gr=' + str(Gr)
    ax.plot(z, U, linestyle='-', color=color_q2d[i], label=label_q2d)
    # Plot validation data
    label_val = 'vetcha ' + 'Gr=' + str(Gr)
    z_val = val_data[val_data.columns[2*i]]
    U_val = val_data[val_data.columns[2*i+1]]
    ax.plot(z_val, U_val, linestyle='--', color=color_val[i], label=label_val)
    # # Get metric of performance:
    # #   root mean squared error over the normalized velocities
    # rmse = np.sqrt(mean_squared_error(y_true=U_val, y_pred=U))
    # print('Root mean squared error of the normalized velocity is', rmse)
    # errors[i+1] = [Ha, Re, Gr, rmse]
    # Update counter
    # Store the simulated case in other directory
    sp.call("mv case case_" + str(Ha) + "_" + str(Re) + "_" + str(Gr), shell=True)
    i += 1

# # Write validation errors file
# pd.DataFrame.from_dict(errors, orient='index',
#                        columns=['Ha', 'Re', 'Gr', 'rmse']).to_csv(
#                        'validation_error_vetcha.csv', index=False)

ax.set_title('Vetcha validation. Ha=50, Re=1e4')
ax.legend(loc='best')
fig.savefig('validation_vetcha.png',format='png',dpi=200)
plt.show()
