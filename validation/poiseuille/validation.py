#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
# import custom modules
sys.path.insert(1, '../../python')
from meshAndGo import meshAndGo
from MHDutils import getLatestTime

# Setup case
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
phys_dict = {
    'rho0'   : 9720,
    'nu'     : 1.54e-7,
    'Cp'     : 189,
    'k'      : 22.36,
    'beta'   : 1.2e-4,
    'sigma'  : 763000
}
mesh_dict = {       # SETUP FOR CYCLIC ONE-CELL CASE
    'Lx'    : 0.1,
    'LxHalf': 0.1/2,
    'Nx'    : 1,
    'c2c_bl'    : 1.05,
    'c2c_bulk'  : 1.001
}
postprocess_dir = 'case/postProcessing/sets/'
postprocess_file = 'line_centreProfile_U.xy'
fig, ax = plt.subplots(figsize=(12,14))

# Run case
meshAndGo(Ha=0, Re=100, Gr=0, Nx=1, Lx=0.1,
    mesh_dict=mesh_dict, tag_dict=tag_dict, phys_dict=phys_dict)

# Get latest time directory name and load simulation data
latest_time = getLatestTime(postprocess_dir)
filename = postprocess_dir + latest_time + '/' + postprocess_file
z, u_q2d, _, _ = np.loadtxt(filename, unpack= True)

# Solution for the Poiseuille Flow between two infinite plates
_, _, _, gradP, _, _ = np.loadtxt('case/output.out', unpack= True)
h = 2*tag_dict['a']
y = np.linspace(0, h)
deltaP = gradP[-1]
mu = phys_dict['rho0']*phys_dict['nu']
u_val = (deltaP*phys_dict['rho0']) /(2*mu) * y * (h-y)

# Plot
ax.plot(y, u_val, label='Poiseuille')
ax.plot(z, u_q2d, linestyle='dotted', color='r', linewidth=3, label='Q2D')
ax.grid(True)
ax.set_ylabel('Velocity (m/s)')
ax.set_xlabel('Channel length (m)')
ax.legend(loc='best')

print("Q2DmhdFOAM's mean velocity:", u_q2d.mean())
print("Poiseuille's mean velocity:", u_val.mean())
fig.savefig('validation_poiseuille.png',format='png',dpi=300)
plt.show()
