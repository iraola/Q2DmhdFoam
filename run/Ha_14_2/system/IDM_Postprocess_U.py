#!/usr/bin/env python3
#import pdb; pdb.set_trace()

#...............................
#Instabilities Detection Method a.k.a IDM
#SCRIPT 2: OpenFOAM Postprocess - Calculates KEC and HV algorithms and returns an instabilities indicator and (optionally) charts.
#AUTHOR: J.Serrat / joaquimserratserra@gmail.com
#SUPERVISION: D.Cambra & E.Iraola
#ENVRIRONMENT: ANT Research Group (UPC)
#...............................

import math
from math import log
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as ml
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import scipy.stats
from decimal import Decimal
import os

#MAIN DATA...................................................................
path = os.getcwd()                    #Current directory path
case = os.path.basename(path)         #Get the current working directory name
field = 'U' #U
fileName = '/postProcessing/probes/0/' + str(field)       #Probed data directory
png = case + 'IDM_chart'               #Setting file names
png2 = case + str(field)               #Field probed
#param = np.loadtxt('MHDinput', unpack=True, skiprows=2)
#Ha = param[0]                         #Hartmann number (optional)
#Re = param[1]                         #Reynolds number (optional)
#Gr = param[2]                         #Grashof number (optional)
U_mean = 2.1e-4                   #nu / l (3D_cylinder)
f1 = 2                                 #Time step calculated                          


#PRAMETERS TO BE CALIBRATED....................................................................
#KEC
maxScale =      4  
powerLimit =    1e4
#HV
param1 =        10
param2 =        10
param3 =        20/.133333333
param4 =        100 

#FLAGS........................................................................
#1D plot
plot_chart =    1                          #0=don't plot / 1=plot & store
#2D plot
plot_chart_2D = 1                       #0=don't plot / 1=plot & store

#PROBED GRID...................................................................
n = 8                               #2**n Probed points per 2**n rows (2**n x n**2 matrix) 
L = 2
Ly= 2                           #Channel width
pp_x = 2.25                           #X coordinate center
pp_y = 1                          #Y coordinate center
pp_z = 1                           #Z coordinate center
delta_x = L/(2**n+1)                    #Distance between probed points
delta_y = Ly/(2**n+1) 
xx = np.arange(-L/2 + pp_x + delta_x/2, L/2 + pp_x - delta_x/2, delta_x)      #Horizontal  points
yy = np.arange(-Ly/2 + pp_y + delta_y/2, Ly/2 + pp_y - delta_y/2, delta_y)                          #Vertical points
X, Y = np.meshgrid(xx, yy)            #Points grid

#READ DATA.................................................................................................
def readProbes(n, path, fileName, f):       #Import data from "probes/0/U" to np.arrays, e.g. a1=readProbes(n, path, fileName1, f1)
     import numpy as np
     text_file = open(path+fileName, 'r')
     lines = text_file.readlines()
     l=lines[2**(n*2)+f]
     x=np.array(l.split()[1:])
     lll=np.array([])
     for k in range(len(x)):
         if (k)%3 == 0:
             lll=np.append(lll , np.array(x[k][1:]).astype(np.float))
     a=np.asmatrix(np.split(lll,2**(n)))
     at=a.transpose()
     att=np.squeeze(np.asarray(at))
     aa=np.squeeze(np.asarray(a))
     return a,aa,at,att
     
#LOG SLOPE.................................................................................................
def log_slope(x1,x2,y1,y2):                 #Calculate log slope from x1,x2,y1,y2
    return (log(y1,10)-log(y2,10))/(log(x2,10)-log(x1,10))

#KOLMOGOROV ENERGY CASCADE ALGORITHM (KEC).................................................................
def KEC(at, L, n, rn, maxScale, powerLimit):                      #at=data array, L=channel width, 2**n probes, rn=row number      
        b = np.abs(at)
        att=np.squeeze(np.asarray(at))
        freq=np.fft.rfftfreq(2**n,L)*2**n
        fft=np.abs(np.fft.rfft(att[rn]/1)/(2**(n-1)))**(2)
        max_fft=np.max(fft)
        index_list=[]
        index_fft=[]
        for p,q in enumerate(fft):
            flag=0
            if q>max_fft/powerLimit and flag==0:
                index_list.append(p)
                index_fft.append(q)                
            else:
                flag=1
        max_list=np.max(np.array(index_list))
        min_list=np.min(np.array(index_list))
        x1=freq[maxScale]
        x2=freq[max_list]
        y1=fft[maxScale]
        y2=fft[max_list]
        slope=-log_slope(x1,x2,y1,y2)
        k1=np.log10(np.array(index_list))
        k2=np.log10(0.5/L*np.array(index_fft))
        slopek, intercept, r_value, p_value, std_err = scipy.stats.linregress(k1[maxScale:], k2[maxScale:])
        return slope,slopek,intercept,r_value**2

#HORITZONAL-VERTICAL ALGORITHM (HV)........................................................................
def HV(field_vector,L,n,slope,interception,r2):             #e.g HV(a1[1],L,n)
    bb = np.fft.rfft2(field_vector) #FFT2
    b = np.abs(bb)
    bbb=b/2**(2*n)
    max_fft2=np.max(bbb)
    l1=[]
    l2=[]
    weak_turb=[]
    strong_turb=[]
    freq=np.fft.rfftfreq(2**n,L)*2**n
    k, m, i = 0, 0, 0
    for k in range(param1):
        if max_fft2<param3*np.max(bbb[k+1][1:]):
            l1.append(k+1)            
    
    if len(l1)!=0:
        for i in range(param1):
            m=0
            if i+1 == l1[m]:
                m+=1
                if np.average(bbb[i+1][1:param2])/max_fft2*param4 > 1:
                    l2.append(i+1)
                
    for p2 in l2:
        strong_turb.append([freq[p2],np.max(bbb[p2][1:])/max_fft2,np.average(bbb[p2][1:param2])/max_fft2])
    l3=list(set(l2)^set(l1))
    for p1 in l3:
        weak_turb.append([freq[p1],np.max(bbb[p1][1:])/max_fft2])
    FFT_file= open("IDM_output_U.txt","w+")
    if weak_turb==[] and strong_turb==[]:
        FFT_file.write('Slope: '+str(slope)+'\n'+'Interception: '+str(interception)+'\n'+'R Squared: '+str(r2)+'\n'+'Fully developed behaviour - No vertical variability. '+'\n'+' For more info plot FFT2'+'\n')
        FFT_file.close() 
    else:
        FFT_file.write(str(slope)+'\n'+'Weak turbulence:' +str(weak_turb)+'\n'+'Strong turbulence:' +str(strong_turb)+'. '+'\n'+'For more info plot FFT2'+'\n'+'DETAILS: Weak turbulence refers to specific periodic eddies [Wavenumber, PEAK/MAX]'+'\n'+'DETAILS: Strong turbulence refers to wide scale range of turbulent structures [Wavenumber, PEAK/MAX, AVG/MAX]')
        FFT_file.close()        

a1=readProbes(n, path, fileName, 2)  #Imported DATA
hh=np.array([])
for ko in range(2**n-1):              #Average 2**n slopes  
    hh0=np.append(hh,KEC(a1[2],Ly,n,ko,maxScale,powerLimit)[1])
    hh1=np.append(hh,KEC(a1[2],Ly,n,ko,maxScale,powerLimit)[2])
    hh2=np.append(hh,KEC(a1[2],Ly,n,ko,maxScale,powerLimit)[3])
HV(a1[1],L,n,np.average(hh0),np.average(hh1),np.average(hh2))          #Execution & storage KEC/HV results


#CHARTS ................................................................................................
if plot_chart == True:
    #APLITUDE LOG SPECTRUM
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, figsize=(9, 9))
    ax0.loglog(np.fft.rfftfreq(2**n,Ly)*2**n, np.abs(np.fft.rfft(a1[3][2**(n-2)])/(2**(n)))/U_mean,'.-')
    ax0.set_title('APLITUDE LOG SPECTRUM')
    ax0.set_ylim([1e-6,1e2])
    ax0.grid()
    ax0.set(xlabel='Wavenumber [$m^{-1}$]', ylabel='Amplitude [Adimentionalized field units]')
    #POWER LOG SPECTRUM
    ax1.loglog(np.fft.rfftfreq(2**n,Ly)*2**n, np.abs(np.fft.rfft(a1[3][2**(n-2)])/(2**(n))/U_mean)**2,'.-')
    ax1.set_title('POWER LOG SPECTRUM')
    ax1.set_ylim([1e-10,1e2])
    ax1.grid()
    ax1.set(xlabel='Wavenumber [$m^{-1}$]', ylabel='Power [Adimentionalized squared field units]')
    #FIELD PROFILE
    ax2.plot(np.arange(-Ly/2+delta_y, Ly/2, delta_y), a1[3][2**(n-2)]/U_mean,'.-')
    ax2.set_title('FIELD PROFILE')
    margin_y = 9/8
    min_chart = np.min(a1[3][2**(n-2)])/U_mean
    max_chart = np.max(a1[3][2**(n-2)])/U_mean
    Amp = max_chart - min_chart
    ax2.set_ylim([min_chart*margin_y,min_chart+Amp*margin_y])
    #ax2.set_ylim([0,0.4])
    ax2.grid()
    ax2.set(xlabel='Channel width [m]', ylabel='Amplitude [Adimentionalized field units]')
    plt.subplots_adjust(left=0.125, bottom=0.075, right=0.9, top=0.975, wspace=0.2, hspace=0.45)
    #plt.show()
    plt.savefig('3D_cylinder_FFT_U',dpi=600)

#2D CHARTS.......................................................................................
if plot_chart_2D == True:
    #Velocity profile
    fig, ax = plt.subplots()
    cs=ax.imshow(a1[3]/U_mean, interpolation='none', cmap='coolwarm',
                   origin='lower', extent=[-Ly/2, Ly/2, -Ly/2, Ly/2])
    ax.set(xlabel='Channel width y-axis [m]', ylabel='Channel width z-axis [m]')
    ax.set_title('FIELD PROFILE [Adimentionalized velocity]')
    cbar = fig.colorbar(cs)
    plt.savefig('3D_cylinder_FieldProfile_U',dpi=600)
    #2D profile spectrum
    fig, ax = plt.subplots(figsize=(6, 3))
    bb = np.fft.fft2(a1[3]) #/U_mean #FFT2
    b = np.abs(bb)
    bbb=b/2**(2*n)
    kl_1=bbb[0:2**(n-1),0:2**(n-1)]
    kl_2=bbb[2**(n-1):2**(n),2**(n-1):2**n]
    kl_3=bbb[0:2**(n-1),2**(n-1):2**n]
    kl_4=bbb[2**(n-1):2**(n),0:2**(n-1)]
    #kl_1+kl_4[::-1]i
    U_mean = np.max(kl_1)
    kl_1 /= U_mean
    freq=np.fft.rfftfreq(2**n,Ly)*2**n
    cs=ax.imshow(kl_1[0:2**(n-5),0:2**(n-4)], interpolation="none", cmap='Spectral',
        origin='lower',extent=[0, freq[2**(n-4)], 0, freq[2**(n-5)]],vmax=np.unique(kl_1)[-2], vmin=0)
    ax.set(xlabel='y-Wavenumber [$m^{-1}$]', ylabel='z-Wavenumber [$m^{-1}$]')
    ax.set_title('FFT2 [Adimentionalized velocity]')
    ticks = np.linspace(0, np.unique(kl_1)[-2], 6)
    cbarf = fig.colorbar(cs, ticks=ticks)
    #cbar.ax.set_yticklabels(['{%.1E}'.format(i) for i in ticks])
    plt.subplots_adjust(left=0.14, bottom=0.2, right=1, top=0.9, wspace=0.0, hspace=0.0)
    #plt.show()
    ax.set_xlim(0, freq[2**(n-4)])
    ax.set_ylim(0, freq[2**(n-5)])
    ax.xaxis.set_minor_locator(AutoMinorLocator(3))
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))
    ax.xaxis.set_major_locator(MultipleLocator(1.5))
    ax.yaxis.set_major_locator(MultipleLocator(1.5))
    # Turn grid on for both major and minor ticks and style minor slightly
    # differently.
    ax.grid(which='major', color='#CCCCCC', linestyle='-')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    
    #ax.set_aspect(1/2, adjustable=None, anchor=None, share=False)

        
    vegetables = kl_1[0:2**(n),0:2**(n)]
    for i1,e1 in enumerate(vegetables):
        for i2,e2 in enumerate(e1):
            if vegetables[i1,i2] > 0.02:
                text = ax.text((i2*1+0.5)/2, (i1*1+0.5)/2, "%1.2f" %vegetables[i1,i2], ha="center", va="center", color="w",fontsize=6)
    
    #text = ax.text(0*1+0.5, 0*1+0.5, 0.0, ha="center", va="center", color="w",fontsize=8)
    ax.annotate("$A_{0,0}$", xy=(0.5, 0.5), xycoords="data",
        xytext=(1.5, 1.5), va="top", ha="center",
        bbox=dict(boxstyle="round", fc="w"),
        arrowprops=dict(arrowstyle="->"))
        
    ax.annotate("$A_{i,0}$", xy=(3,0.25 ), xycoords="data",
        xytext=(101.33/6.666/2, 0.38), va="top", ha="center",
        bbox=dict(boxstyle="round", fc="w"),
        arrowprops=dict(arrowstyle="<|-|>"))
                      
    ax.annotate("$A_{0,j}$", xy=(0.25, 1), xycoords="data",
        xytext=(0.25, 52/6.666/2), va="top", ha="center",
        bbox=dict(boxstyle="round", fc="w"),
        arrowprops=dict(arrowstyle="<|-|>"))
        
    plt.savefig('3D_cylinder_FFT2_U',dpi=600)
                    
