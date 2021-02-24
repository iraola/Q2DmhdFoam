#!/usr/bin/env python3
#import pdb; pdb.set_trace()
import MHDutils as utils
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
import os

######################
### INITIALIZATIONS
######################
# Get main parameters from input file
param = np.loadtxt('MHDinput', unpack=True, skiprows=2)
Ha = param[0]
Re = param[1]
Gr = param[2]

# Define main tag dictionary
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
    'g'   : -9.81
}

# Update main dictionary with physical properties
phys_dict = {
    'rho0'   : 9720,
    'nu'     : 1.54e-7,
    'Cp'     : 189,
    'k'      : 22.36,
    'beta'   : 1.2e-4,
    'sigma'  : 763000
}
tag_dict.update(phys_dict)

# Define directory names
base_dir = 'baseCaseTagawa'
curr_dir = 'tagawa'
# Create new case directory
sp.call('cp -r ' + base_dir + ' ' + curr_dir, shell=True)

######################
### CASE CONDITIONS
######################
# Calculations from the parametrization of the problem and assign
# them to the dictionaries
B = Ha / (2*tag_dict['b']) * np.sqrt(tag_dict['nu'] * tag_dict['rho0'] /
    tag_dict['sigma'])
Ux = Re * tag_dict['nu'] / tag_dict['a']
# Calculate delta_T from the input Grashof number
delta_T = Gr * tag_dict['nu']**2 / (abs(tag_dict['g']) * tag_dict['beta']
    * (2*tag_dict['a'])**3 )
# Get hot and cold temperature for walls
Th = delta_T/2
Tc = -delta_T/2
# Tag and print values
tag_dict['B'] = B
tag_dict['Ux'] = Ux
tag_dict['Th'] = Th
tag_dict['Tc'] = Tc
print('Magnetic field B = {} for Ha = {}'.format(tag_dict['B'],Ha))
print('Mean velocity U = {} for Re = {}'.format(tag_dict['Ux'],Re))
print('Temperature difference is {} for Grashof = {}'.format(delta_T,Gr))

######################
### MESH
######################

l_side = tag_dict['a']/np.sqrt(Ha)
N_bl = 25
G_bl, _, cmax = utils.lenC2CN(l_side, N_bl)

print('Ha = {}'.format(Ha))
print('SIDE BOUNDARY LAYER')
print('side boundary length = {}'.format(l_side))
print('G = {}'.format(G_bl))
print('cmax = {}'.format(cmax))

# This is HALF the length of the bulk
l_bulk = tag_dict['a'] - l_side
N_bulk, G_bulk, _ = utils.lenCminC2C(l_bulk, cmax)
print('HALF BULK')
print('bulk half length = {}'.format(l_bulk))
print('N = {}'.format(N_bulk))
print('G = {}'.format(G_bulk))

mesh_dict = {
    'Lx'        : 200*tag_dict['a'],
    'LxHalf'    : 200*tag_dict['a'] / 2,
    'LyBulk'    : tag_dict['a'],
    'LyBL'      : tag_dict['a']-l_side,
    'LyNegBulk' : -tag_dict['a'],
    'LyNegBL'   : -(tag_dict['a']-l_side),
    'Nx'        : 100,
    'NyBulk'    : N_bulk,
    'NyBL'      : N_bl,
    'GyBulk'    : G_bulk,
    'GyBulkInv' : 1/G_bulk,
    'GyBL'      : G_bl,
    'GyBLinv'   : 1/G_bl
}
# Combine the mesh dictionary with the previous dictionary
tag_dict.update(mesh_dict)
# Replace all tags
utils.replaceTags(tag_dict, curr_dir)


######################
### RUN CASE
######################
os.chdir(curr_dir)
sp.call("blockMesh > log.blockMesh", shell=True)
sp.call("Q2DmhdFoam > log.Q2DmhdFoam", shell=True)
sp.call("sample > log.sample", shell=True)

######################
### PLOT & VALIDATION
######################
# Analytical solution
T = np.linspace(Th,Tc)
theta = (T - 1/2*(Th + Tc))/(Th-Tc)
y_analytic = np.linspace(-tag_dict['a'],tag_dict['a'])
W_analytic = Gr/(2*Ha) * theta
fig, ax = plt.subplots(figsize=(12,12))
ax.plot(y_analytic[10:-10]/(2*tag_dict['a']), W_analytic[10:-10],
    label='Bulk analytical solution')

# Q2DmhdFoam solution
postProcess_dir = 'postProcessing/sets/'
latestTime = utils.getLatestTime(postProcess_dir)
raw_data_file = postProcess_dir + latestTime + '/line_centreProfile_U.xy'
y, Ux, Uy, Uz = np.loadtxt(raw_data_file, unpack=True)
y = y - tag_dict['a']
W = Ux / (tag_dict['nu'] / (2*tag_dict['a']))
ax.plot(y/(2*tag_dict['a']), W, label='Q2DmhdFoam')
# Print max velocity
print('The Q2D maximum velocity is : ', np.max(W))

# Plot derivative of the Q2D solution
# Note: Use dimensionless Y. Y = y / l
y_plot = np.zeros(len(y)-1)
Y = y / (2*tag_dict['a'])
dW = np.zeros(len(W)-1)
for i in range(len(W) - 1):
    dW[i] = (W[i+1] - W[i]) / (Y[i+1] - Y[i])
    y_plot[i] = (Y[i] + Y[i+1])/2
fig1, ax1 = plt.subplots(figsize=(12,12))
ax1.plot(y_plot, -dW, label='Q2D dW/dY')
ax1.set_ylim(bottom=0)

# Write out tagawa values
out_file = '../tagawa.csv'
out_line = str(Ha) + ' ' + str(np.max(W)) + ' ' + str(
    -dW[int((len(dW)-1)/2)]) + '\n'
# Check if the file exists
if not os.path.exists(out_file):
    text = out_line
else:
    with open(out_file, 'r') as file:
        text = file.read()
    text = text + out_line
# Write file
with open(out_file, 'w') as file:
    file.write(text)

# Plot config
ax.set_xlabel('Dimensionless length, $X=x/l$')
ax.set_ylabel('Dimensionless velocity, $W$')
ax.legend(loc='upper right')
plt.show()
fig.savefig('tagawaVal.eps',format='eps',dpi=1000)
