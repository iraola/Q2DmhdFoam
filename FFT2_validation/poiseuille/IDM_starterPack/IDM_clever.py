import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#import matplotlib.animation as animation
import sys
import os
#import pandas as pd
import subprocess
from decimal import Decimal



def returnValue(file_, find_str):
#    find_str = "pressure gradient"                    # String to find
    with open(file_, "r") as f:
        f.seek (0, 2)           # Seek @ EOF
        fsize = f.tell()        # Get Size
        f.seek (max (fsize-1024, 0), 0) # Set pos @ last n chars
        lines = f.readlines()       # Read to end
        lines = lines[-200:]    # Get last 10 lines
        for line in lines:
            if find_str in line:
                value = line.rpartition(find_str)[2]
                value = value.rstrip()
                value = value.rstrip(".")
                value = value.rstrip("J/s")
                break
    return value


wd = os.getcwd()
directoris = os.listdir(wd)
print ("Current folder is " + str(wd))
###########################################################
##### DEFINICION CASOS:

Ha = []
Re = []
Gr = []
N  = []

IDM_value=-3
slope_vort = []
slope_U = []
Ca1_ = []
Ca2_ = []
Ca3_ = []

# Entramos en las carpetas de los casos y leemos IDM_output.txt:
for carpeta in directoris:
    if "case_" in carpeta:
        Ca1_ = carpeta.split("Re")[1:]
        Ca2_ = Ca1_[0].split("Gr")[1:]
        Ca3_ = Ca2_[0].split("Ha")[1:]
        Re_ = float(Ca1_[0].split("Gr")[0])
        Gr_ = float(Ca2_[0].split("Ha")[0])
        Ha_ = float(Ca3_[0].split("Ha")[0])
        N_  = Ha_**2/Re_
        file_interest_vort   = carpeta + '/IDM_output.txt'
        file_interest_U   = carpeta + '/IDM_output_U.txt'
        file_log = str("case_"+str(int(Ha_))+"_"+str(int(Re_))+"_"+str(float(Gr_))+'/log')
        try:
            first_line_vort = subprocess.check_output(["sed","-n","1p", file_interest_vort])
            first_line_U = subprocess.check_output(["sed","-n","1p", file_interest_U])
            slope_vort_ = float(first_line_vort.split()[0])
            slope_U_ = float(first_line_U.split()[0])
            Ha.append(Ha_)
            Re.append(Re_)
            Gr.append(Gr_)
            N.append(N_)
            slope_vort.append(slope_vort_)
            slope_U.append(slope_U_)
        except:
            print('First except '+carpeta)
# ----------------------------------------------------------------------

print('Re \t Gr \t \t Ha \t N \t \t slope_vort \t slope_U  ')
print('===============================================================================')
for i in range(len(Ha)):
    print(str(int(Re[i]))+' \t '+'%.1E' % Decimal(Gr[i])+' \t '+str(int(Re[i]))+' \t '+'%.1E' % Decimal(N[i])+' \t '+"%1.4f" %slope_vort[i]+' \t '+"%1.4f" %slope_U[i] )
print("")

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
#ax = plt.axes(projection='3d')
instability='^'
stability='o'
"""
# Add x, y gridlines
ax.grid(b = True, color ='grey',
        linestyle ='-.', linewidth = 0.3,
        alpha = 0.2)
 """

# Creating color map
for k in range(len(Ha)):
    xs = Ha[k]
    ys = Gr[k]
    zs = Re[k]
    if slope_vort[k]<IDM_value:
        ax.scatter(xs, ys, zs,c='b', marker=stability)
        #ax.scatter(xs, ys, zs,c=abs(slope_vort[k]), cmap='viridis', linewidth=0.5 , marker=stability)
    else:
        ax.scatter(xs, ys, zs,c='r', marker=instability)
        #ax.scatter(xs, ys, zs,c=abs(slope_vort[k]), cmap='viridis', linewidth=0.5 , marker=instability)
  
plt.title("simple 3D scatter plot")     
ax.set_xlabel('Ha')
ax.set_ylabel('Gr')
ax.set_zlabel('Re')
plt.show()
