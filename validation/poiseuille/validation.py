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
postprocess_dir = 'case/postProcessing/sets/'
postprocess_file = 'line_centreProfile_U.xy'
fig, ax = plt.subplots(figsize=(12,6))

# Run case
meshAndGo(Ha=0, Re=100, Gr=0, Nx=1, Lx=0.1,
    tag_dict=tag_dict, phys_dict=phys_dict)

# Get latest time directory name and load simulation data
latest_time = getLatestTime(postprocess_dir)
filename = postprocess_dir + latest_time + '/' + postprocess_file
z, U, _, _ = np.loadtxt(filename, unpack= True)
ax.plot(z, U, label='Q2D')

# Solution for the Poiseuille Flow between two infinite plates
_, _, _, gradP, _, _ = np.loadtxt('case/output.out', unpack= True)
h = 2*tag_dict['a']
y = np.linspace(0, h)
deltaP = gradP[-1]
mu = phys_dict['rho0']*phys_dict['nu']
u = deltaP/(2*mu) * y * (h-y)
ax.plot(y,u, label='Poiseuille')
ax.legend(loc='best')
plt.show()
