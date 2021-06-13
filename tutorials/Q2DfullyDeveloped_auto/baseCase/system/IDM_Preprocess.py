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
n = /n/                               #2**n Probed points per 2**n rows (2**n x n**2 matrix) 
L = /a/*2                             #Channel width
pp_x = /x0/                           #X coordinate center
pp_y = /y0/                           #Y coordinate center
delta = L/(2**n+1)                    #Distance between probed points
xx = np.arange(-L/2 + pp_x + delta/2, L/2 + pp_x - delta/2, delta)      #Horizontal  points
yy = np.arange(-L/2 + pp_y + delta/2, L/2 + pp_y - delta/2, delta)                          #Vertical points
X, Y = np.meshgrid(xx, yy)            #Points grid


def coordinates(x,y):                 #Convert probed points into STR data e.g. (0, -0.015820312500000558, 0.017578124999999126)
    l = []
    for i_y in range(len(y)):
        for i_x in range(len(x)):
            l += ['('+str(x[i_x])+' '+str(y[i_y])+' '+str(0)+')']
    return l

fileName = 'points'
with open(fileName, 'w') as f1:        #Store points to "fileName" file
    for item in coordinates(xx,yy):
        f1.write("%s\n" % item)
f1.close()        
#END...........................................................................
#NEXT STEP: cat controlDict1 "fileName" controlDict2 >> controlDict