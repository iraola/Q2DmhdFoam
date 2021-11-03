#!/usr/bin/env python3
import numpy as np
from numpy import pi, sin, sinh
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
    'Nx'    : 1}
# Get physical properties used later
a = tag_dict['a']
b = tag_dict['b']
rho = phys_dict['rho0']
nu = phys_dict['nu']
mu = rho * nu
# Procedural variables
tol = 1e-15
postprocess_dir = 'case/postProcessing/sets/'
postprocess_file = 'line_centreProfile_U.xy'
fig, ax = plt.subplots(figsize=(12,6))
fig3d, ax3d = plt.subplots(figsize=(12,6))
ax3d = plt.axes(projection='3d')

# Run case
meshAndGo(Ha=0, Re=100, Gr=0, Nx=1, Lx=0.1,
    mesh_dict=mesh_dict, tag_dict=tag_dict, phys_dict=phys_dict)

# Get latest time directory name and load simulation data
latest_time = getLatestTime(postprocess_dir)
filename = postprocess_dir + latest_time + '/' + postprocess_file
z, u_q2d, _, _ = np.loadtxt(filename, unpack= True)
ax.plot(z, u_q2d, label='Q2D')

# Solution for the Poiseuille Flow in a square pipe. It needs recursion
# # Important: remember that pressure in OF is divided by density
n_grid = 100
_, _, _, gradP, _, _ = np.loadtxt('case/output.out', unpack= True)
deltaP = gradP[-1]
G = deltaP * rho
h = 2 * a
l = 2 * b
y = np.linspace(0, h, n_grid)
z = np.linspace(0, l, n_grid)
yv, zv = np.meshgrid(y, z)
u_val = G / (2*mu) * yv * (h - yv)
u_sum_temp = np.zeros((n_grid, n_grid))
n = 1
end_flag = False
mod_old = 0.1

while not end_flag:
    print('\nfor n =', n)
    beta = (2*n - 1) * pi / h
    print('beta:', beta)
    u_sum_temp += 1 / (2*n - 1)**3 * (sinh(beta * zv) + sinh(beta * (l-zv))) \
                                    / sinh(beta * l) * sin(beta * yv)
    mod_new = u_sum_temp.mean().mean()
    print('mod_new:', mod_new)
    print('mod_old:', mod_old)
    if abs((mod_new - mod_old)/mod_old) < tol:
        end_flag = True
    mod_old = mod_new.copy()
    n += 1

ax.plot(y, u_val[0,:], label='Poiseuille 1D')
u_val = u_val - 4*G*h**2 / (mu*pi**3) * u_sum_temp
ax.plot(y, u_val[int(u_val.shape[0]/2), :], label='Poiseuille 2D')
ax3d.plot_surface(yv, zv, u_val)
ax.legend(loc='best')

print("Q2DmhdFOAM's mean velocity:", u_q2d.mean())
print("Poiseuille's mean velocity:", u_val.mean())
fig.savefig('validation_poiseuille.png',format='png',dpi=200)
plt.show()
