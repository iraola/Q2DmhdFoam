#!/usr/bin/env python3
#import pdb; pdb.set_trace()

#...............................
#Instabilities Detection Method a.k.a IDM
#SCRIPT 1: OpenFOAM Preprocess - Calculates the probed points annd modifies the controlDict for probes funcionality.
#AUTHOR: J.Serrat / joaquimserratserra@gmail.com
#SUPERVISION: D.Cambra & E.Iraola
#ENVRIRONMENT: ANT Research Group (UPC)
#...............................

import numpy as np
import os

#PROBED GRID...................................................................
n = 8                               #2**n Probed points per 2**n rows (2**n x n**2 matrix) 
L = 0.1
Ly=0.15                           #Channel width
pp_x = 0.05                           #X coordinate center
pp_y = 0                           #Y coordinate center
pp_z = 0.05                           #Z coordinate center
delta_x = L/(2**n+1)                    #Distance between probed points
delta_y = Ly/(2**n+1) 
xx = np.arange(-L/2 + pp_x + delta_x/2, L/2 + pp_x - delta_x/2, delta_x)      #Horizontal  points
yy = np.arange(-Ly/2 + pp_y + delta_y/2, Ly/2 + pp_y - delta_y/2, delta_y)                          #Vertical points
X, Y = np.meshgrid(xx, yy)            #Points grid


def coordinates(x,y):                 #Convert probed points into STR data e.g. (0, -0.015820312500000558, 0.017578124999999126)
    l = []
    for i_y in range(len(y)):
        for i_x in range(len(x)):
            l += ['('+str(x[i_x])+' '+str(y[i_y])+' '+str(0.05)+')']
    return l

fileName = 'points'
with open(fileName, 'w') as f1:        #Store points to "fileName" file
    for item in coordinates(xx,yy):
        f1.write("%s\n" % item)
f1.close()        
#END...........................................................................
#NEXT STEP: cat controlDict1 "fileName" controlDict2 >> controlDict