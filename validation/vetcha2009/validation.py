#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
# import custom modules
sys.path.insert(1, '../../python')
from meshAndGo import meshAndGo

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
mesh_dict = {
    'Lx'    : 0.1,
    'LxHalf': 0.1/2,
    'Nx'    : 1
}
# Run case
meshAndGo(Ha=100, Re=100, Gr=1e4,
    mesh_dict=mesh_dict, tag_dict=tag_dict, phys_dict=phys_dict)

# Solution for the Poiseuille Flow between two infinite plates
h = 2*tag_dict['a']
y = np.linspace(0, h)
deltaP = 1
mu = phys_dict['rho0']*phys_dict['nu']
u = deltaP/(2*mu) * y * (h-y)
plt.plot(y,u)
plt.show()
