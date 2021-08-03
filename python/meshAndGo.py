#!/usr/bin/env python3
#import pdb; pdb.set_trace()
import MHDutils as utils
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
import os

def meshAndGo(Ha, Re, Gr, volumetric_heat=True, write=False, plot=False,
                verbose=False, **kwargs):
    '''
    Setup, mesh, run and plot results for a Q2DmhdFOAM case

    Parameters
    ----------

    '''
    print('\n\nStarting simulation for Ha={}, Re={}, Gr={}\n'
        .format(Ha, Re, Gr))
    ######################
    ### INITIALIZATIONS
    ######################
    # Get main parameters from input file
    # param = np.loadtxt('MHDinput', unpack=True, skiprows=2)
    # Ha = param[0]
    # Re = param[1]
    # Gr = param[2]

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
        'Th'  : 0,
        'Tc'  : 0,
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
    # Optional arguments for changing default dicts using kwargs
    #   NOTE: if the user provides parameters that already exist in the dicts,
    #   they will be overwritten. The rest will remain the same.
    #   This is very useful if you want to use all defaults except a couple of
    #   variables
    if 'tag_dict' in kwargs:
        tag_dict.update(kwargs['tag_dict'])
    if 'phys_dict' in kwargs:
        phys_dict.update(kwargs['phys_dict'])

    # Add all tags into same dict
    tag_dict.update(phys_dict)

    # Get dict variables for easy handling
    a = tag_dict['a'];      b = tag_dict['b']
    nu = tag_dict['nu'];    rho = tag_dict['rho0'];     sigma = tag_dict['sigma']
    g = tag_dict['g'];      beta = tag_dict['beta'];    k = tag_dict['k']
    m = tag_dict['m']

    # Define directory names
    base_dir = 'baseCase'
    curr_dir = 'case'
    # Create new case directory and remove old if exists
    if os.path.exists(curr_dir):
        sp.call('rm -r ' + curr_dir, shell=True)
    sp.call('cp -r ' + base_dir + ' ' + curr_dir, shell=True)

    ######################
    ### CASE CONDITIONS
    ######################
    # Calculations from the parametrization of the problem and assign
    # them to the dictionaries
    # FROM HARTMANN:
    B = Ha / (2 * a) * np.sqrt(nu * rho / sigma)

    # FROM REYNOLDS:
    Ux = Re * nu / a

    # FROM GRASHOF:
    delta_T = Gr * nu**2 / (abs(g) * beta * (2 * a)**3)
    print('delta_T', delta_T)
    if volumetric_heat:
        # Calculate exponential deposition q0 (tepot case)
        # - Average heat deposition in W/m3 using Fourier's law
        S0 = delta_T * k / (a**2)
        print('S0', S0)
        # - Calc. q0 for q = q0 * exp(-m * y) as in iniSourceT
        #   Using S0 = 1 / (2*a) * integrate(q from -a to a)
        q0 = 2 * m * S0 / (1 - np.exp(-2 * m))
        # - Do not divide by rho*cp to go from W/m3 to K/s since Q2DmhdFOAM
        #   already understands sourceT comes in W/m3 and divides it by
        #   rho*cp inside the solver
        tag_dict['q0'] = q0
    else:
        # Get hot and cold temperature for walls
        Th = delta_T/2
        Tc = -delta_T/2
        tag_dict['Th'] = Th
        tag_dict['Tc'] = Tc

    # Tag and print values
    tag_dict['B'] = B
    tag_dict['Ux'] = Ux
    print('Magnetic field is B = {} for Ha = {}'.format(B, Ha))
    print('Mean velocity is U = {} for Re = {}'.format(Ux, Re))
    print('Temperature difference is {} for Grashof = {}'.format(delta_T, Gr))

    ######################
    ### MESH
    ######################

    # Default constants (can be ovewritten by user's mesh_dict)
    Nx = 100    # Nx : Number of cells in the direction of the flow
    Lx = 20     # Lx : Proportion of the length in the direction of the flow
                # with respect to the half-length 'a' of the channel
    # Check if user provides cell-to-cell ratios for meshing
    if 'mesh_dict' in kwargs:
        try:
            c2c_bl = kwargs['mesh_dict']['c2c_bl']
            c2c_bulk = kwargs['mesh_dict']['c2c_bulk']
        except:
            print('Please provide "c2c_bl" and "c2c_bulk" in "mesh_dict"')
    else:
        c2c_bl = 1.05
        c2c_bulk = 1.1

    # Mesh MHD calculations
    if Ha == 0:
        l_side = 3 * a / 10
    else:
        l_side = 3 * a / np.sqrt(Ha)
    N_bl = 25
    G_bl, _, cmax = utils.lenC2CN(l_side, N_bl, c2c_bl)
    # This is HALF the length of the bulk
    l_bulk = a - l_side
    N_bulk, G_bulk, _ = utils.lenCminC2C(l_bulk, cmax, c2c_bulk)
    if verbose:
        print('Ha = {}'.format(Ha))
        print('SIDE BOUNDARY LAYER')
        print('side boundary length = {}'.format(l_side))
        print('G = {}'.format(G_bl))
        print('cmax = {}'.format(cmax))
        print('HALF BULK')
        print('bulk half length = {}'.format(l_bulk))
        print('N = {}'.format(N_bulk))
        print('G = {}'.format(G_bulk))

    mesh_dict = {
        'Lx'        : Lx*a,
        'LxHalf'    : Lx*a / 2,
        'LyBulk'    : a,
        'LyBL'      : a-l_side,
        'LyNegBulk' : -a,
        'LyNegBL'   : -(a-l_side),
        'Nx'        : Nx,
        'NyBulk'    : N_bulk,
        'NyBL'      : N_bl,
        'GyBulk'    : G_bulk,
        'GyBulkInv' : 1/G_bulk,
        'GyBL'      : G_bl,
        'GyBLinv'   : 1/G_bl
    }

    # Again, overwrite values if the user provides them:
    if 'mesh_dict' in kwargs:
        mesh_dict.update(kwargs['mesh_dict'])

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

    if plot and not volumetric_heat:
        T = np.linspace(Th,Tc)
        theta = (T - 1/2 * (Th + Tc)) / (Th - Tc)
        y_analytic = np.linspace(-a,a)
        W_analytic = Gr / (2 * Ha) * theta

        # Q2DmhdFoam solution
        postProcess_dir = 'postProcessing/sets/'
        latestTime = utils.getLatestTime(postProcess_dir)
        raw_data_file = postProcess_dir + latestTime + '/line_centreProfile_U.xy'
        y, Ux, Uy, Uz = np.loadtxt(raw_data_file, unpack=True)
        y = y - a
        W = Ux / (nu / (2*a))
        # Print max velocity
        print('The Q2D maximum velocity is : ', np.max(W))

        # Plot derivative of the Q2D solution
        # Note: Use dimensionless Y. Y = y / l
        y_plot = np.zeros(len(y)-1)
        Y = y / (2*a)
        dW = np.zeros(len(W)-1)
        for i in range(len(W) - 1):
            dW[i] = (W[i+1] - W[i]) / (Y[i+1] - Y[i])
            y_plot[i] = (Y[i] + Y[i+1])/2

        # PLOT
        fig, ax = plt.subplots(figsize=(12,12))
        ax.plot(y_analytic[10:-10]/(2*a), W_analytic[10:-10],
            label='Bulk analytical solution')
        ax.plot(y/(2*a), W, label='Q2DmhdFoam')
        fig1, ax1 = plt.subplots(figsize=(12,12))
        ax1.plot(y_plot, -dW, label='Q2D dW/dY')
        ax1.set_ylim(bottom=0)
        # Plot config
        ax.set_xlabel('Dimensionless length, $X=x/l$')
        ax.set_ylabel('Dimensionless velocity, $W$')
        ax.legend(loc='upper right')
        fig.savefig('results.eps',format='eps',dpi=1000)
        plt.show()

    if write:
        # Write out tagawa values
        out_file = '../case.csv'
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

    # Return to original directory
    os.chdir('..')
