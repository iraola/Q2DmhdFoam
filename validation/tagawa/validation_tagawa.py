import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
import pandas as pd
import sys
# import custom modules
sys.path.insert(1, '../../python')
from meshAndGo import meshAndGo
from MHDutils import getLatestTime

def plot_csv(file='tagawa.csv'):
    Ha, maxU, gradU = np.loadtxt(file, unpack=True)
    maxU = np.array(maxU)
    gradU = np.array(gradU)
    Ha = ['Ha='+str(int(ha)) for ha in Ha]
    true_maxU = np.array([16.43, 4.073, 1.123])
    true_gradU = np.array([50.00, 10.0, 2.5])
    rel_error_grad = abs((gradU-true_gradU)/true_gradU)*100
    rel_error_maxU = abs((maxU-true_maxU)/true_maxU)*100

    fig_gradU, ax_gradU = plt.subplots(figsize=(12,12))
    X = np.arange(3) # for ax.bar() purposes
    ax_gradU.bar(X + 0.0, gradU, color = 'b', width = 0.25, label='Analytical Gradient')
    ax_gradU.bar(X + 0.25, true_gradU, color = 'g', width = 0.25, label='Q2D Gradient')
    ax_gradU.bar(X + 0.5, rel_error_grad, color = 'r', width = 0.25, label='relative error')
    ax_gradU.set_xticks(X)
    ax_gradU.set_xticklabels(Ha)
    ax_gradU.set_title('Maximum Velocity')

    fig_maxU, ax_maxU = plt.subplots(figsize=(12,12))
    X = np.arange(3) # for ax.bar() purposes
    ax_maxU.bar(X + 0.0, maxU, color = 'b', width = 0.25, label='tagawa max U')
    ax_maxU.bar(X + 0.25, true_maxU, color = 'g', width = 0.25, label='Q2D max U')
    ax_maxU.bar(X + 0.5, rel_error_grad, color = 'r', width = 0.25, label='relative error')
    ax_maxU.set_xticks(X)
    ax_maxU.set_xticklabels(Ha)
    ax_maxU.set_title('Maximum Velocity')

    ax_gradU.grid()
    ax_gradU.legend(loc='upper right')
    ax_maxU.grid()
    ax_maxU.legend(loc='upper right')
    #ax.bar(X + 0.25, data[1], color = 'g', width = 0.25)
    #ax.bar(X + 0.50, data[2], color = 'r', width = 0.25)

    plt.show()

def plot_bars(values1, values2, descript, labels, values3=None):
    '''
    Plots bars and relative errors with comparation purposes
    '''
    values1 = np.array(values1)
    values2 = np.array(values2)
    rel_error = abs((values1-values2)/values1)*100

    fig, ax = plt.subplots(figsize=(12,12))
    X = np.arange(3) # for ax.bar() purposes
    ax.bar(X - 0.25, values1, color = 'b', width = 0.25, label=labels[0])
    ax.bar(X + 0.0, values2, color = 'g', width = 0.25, label=labels[1])
    ax.bar(X + 0.25, rel_error, color = 'r', width = 0.25, label='relative error')
    ax.set_xticks(X)
    ax.set_xticklabels(descript)
    ax.grid()
    ax.legend(loc='upper right')
    plt.show()

def plot_plot(values1, values2, x, labels, **kwargs):
    '''
    Plots values and relative errors with comparation purposes
    '''
    values1 = np.array(values1)
    values2 = np.array(values2)
    rel_error = abs((values1-values2)/values1)*100

    fig, ax = plt.subplots(figsize=(12,12))
    ax.plot(x, values1, color = 'b', label=labels[0])
    ax.plot(x, values2, color = 'g', label=labels[1])
    ax.plot(x, rel_error, color = 'r', label='relative error (%)')
    ax.grid()
    ax.legend(loc='best')
    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'])
    if 'ylim_min' in kwargs:
        ax.set_ylim(ymin=kwargs['ylim_min'])
    fig.savefig('validation_tagawa.png',format='png',dpi=1000)
    plt.show()

def validation_tepot_AR(Ha, Re, Gr, AR_list):
    '''
    '''
    # Needed tag dicts (because we'll modify them)
    tag_dict = {'B'   : '?',
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
    phys_dict = {'rho0'   : 9720,
        'nu'     : 1.54e-7,
        'Cp'     : 189,
        'k'      : 15.14, # This one 'k' is changed due to discrepancies
        'beta'   : 1.2e-4,
        'sigma'  : 763000}
    maxW = np.zeros(len(AR_list))
    maxW_th = maxW.copy()
    i = 0
    # LOOP
    sp.call("rm -r -f tagawa_*", shell=True)
    for AR in AR_list:
        tag_dict['b'] = tag_dict['a'] * AR
        maxW[i], maxW_th[i] = meshAndGo(Ha, Re, Gr,
                         tag_dict=tag_dict, phys_dict=phys_dict)
        print('\n' + str(maxW[i]) + '\n' )
        i = i + 1
        # Save directory for analysis
        sp.call("mv tagawa tagawa_" + str(AR), shell=True)

    # READ TEPOT VALUES
    _, maxU_tepot = np.loadtxt('tepot_tagawa_100.out', skiprows=1, unpack=True)
    maxW_tepot = maxU_tepot / (tag_dict['nu'] / (2*tag_dict['a']))
    # PLOT
    plot_plot(maxW_tepot, maxW, x=AR_list,
        labels=['tepot max. velocity', 'Q2D max. velocity'], xlabel='AR',
        ylim_min=0)

#AR_list = [1.5, 3, 20./3., 10, 15, 20] # 20/3 = 6.6666...

### INITIALIZATIONS
validation_file = 'tagawa_tab'
postprocess_dir = 'case/postProcessing/sets/'
postprocess_file = 'line_centreProfile_U.xy'
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
    'k'      : 22.36,
    'beta'   : 1.2e-4,
    'sigma'  : 763000}
mesh_dict = {
    'Lx'    : 30,   # x200 times the '2a' length = 200*0.075*2
    'LxHalf': 30/2,
    'Nx'    : 100,
    'c2c_bl'    : 1.05,
    'c2c_bulk'  : 1.001}
# Read validation data
Re = 0
Gr = 1e4
data_val = pd.read_csv(validation_file, engine='python', delimiter='\s',
        names=['Ha', 'maxU', 'Ugrad'], header=None, skiprows=2)
Ha_list = data_val['Ha'].to_numpy()
maxW_th = data_val['maxU'].to_numpy()


### LOOP
# Remove old directories present
sp.call("rm -r -f case_*", shell=True)
fig1, ax1 = plt.subplots(figsize=(12,6))
maxW_q2d = np.zeros(len(Ha_list))
i = 0
for Ha in Ha_list:
    # CAREFUL! USE Ha/2 BECAUSE OF TAGAWA'S DEFINITION OF Ha DIFFERENT THAN THE
    # STANDARD
    meshAndGo(Ha/2, Re, Gr, volumetric_heat=False,
        mesh_dict=mesh_dict, tag_dict=tag_dict, phys_dict=phys_dict)
    # Get latest time directory name and load simulation data
    latest_time = getLatestTime(postprocess_dir)
    filename = postprocess_dir + latest_time + '/' + postprocess_file
    # Plot sim. data
    z, U, _, _ = np.loadtxt(filename, unpack= True)
    # Get dimensionless velocity
    W = U / phys_dict['nu'] * 2 * tag_dict['a']
    z = (z - tag_dict['a']) / tag_dict['a']
    # Store maximum dimensionless velocity
    maxW_q2d[i] = np.max(W)
    print('\nMaximum dimensionless velocity: ' + str(maxW_q2d[i]) + '\n' )
    # Plot simulation profile
    ax1.plot(z, W, label='Q2D Ha='+str(Ha))
    # Store the simulated case in other directory
    sp.call("mv case case_" + str(Ha), shell=True)
    i += 1
ax1.legend(loc='best')
fig1.savefig('validation_tagawa_profiles.png', format='png', dpi=300)

# PLOT
plot_plot(maxW_th, maxW_q2d, x=Ha_list,
    labels=['validation max. velocity', 'Q2D max. velocity'], xlabel='Ha',
    ylim_min=0)
plt.show()
