#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 12:38:17 2017

@author: max
"""
import re
import numpy as np
import os
import math
from math import ceil,e,pi
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from scipy.linalg import solve
from subprocess import call
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)



####################################################################
####################################################################
## SPECIFY DENSITY AND TEMPERATURE PROFILES
####################################################################
####################################################################
ITER=1
if ITER==0:
    ni0 = 3.629E19
    ni9 = 1.523E19
    ni_sep = 0.3E19
    ni_halo = 1E16
    ni_dp = 1E17
    nu_ni = 3.0
    
    ne0 = 3.629E19
    ne9 = 1.523E19
    ne_sep = 0.3E19
    ne_halo = 1E16
    ne_dp = 1E17
    nu_ne = 2.5
    
    Ti0 = 8.54
    Ti9 = 1.56
    Ti_sep = 0.3
    Ti_halo = 0.003
    Ti_dp = 0.03
    nu_Ti = 3.5
    
    Te0 = 3.394
    Te9 = 1.158
    Te_sep = 0.06
    Te_halo = 0.006
    Te_dp = 0.01
    nu_Te = 2.5
    
    pedloc = 0.9
    
    R0_a = 1.67
    kappa = 1.5
    
    impfrac = 0.03 # nz = impfrac * ni
    machine = 'DIII-D'
else:
    ni0 = 0.4E20
    ni9 = 0.35E20
    ni_sep = 0 #0.2E20
    ni_halo = 1E16
    ni_dp = 1E17
    nu_ni = 3.0
    
    ne0 = 1E20
    ne9 = 0.7E20
    ne_sep = 0 #0.4E20
    ne_halo = 1E16
    ne_dp = 1E17
    nu_ne = 2.5
    
    Ti0 = 23*1.2
    Ti9 = 4*1.2
    Ti_sep = 0
    Ti_halo = 0.003
    Ti_dp = 0.03
    nu_Ti = 3.5
    
    Te0 = 24*1.2
    Te9 = 4*1.2
    Te_sep = 0
    Te_halo = 0.006
    Te_dp = 0.01
    nu_Te = 2.5
    
    pedloc = 0.9
    
    R0_a = 1.67
    kappa = 1.5
    
    impfrac = 0.1 # nz = impfrac * ni
    machine = 'ITER'

numpts = 200
rhovals = np.linspace(0,1,numpts)

ni     = np.where(rhovals<0.9,(ni0-ni9)*(1-rhovals**2)**nu_ni + ni9,(ni_sep-ni9)/(0.1)*(rhovals-0.9)+ni9)
Ti     = np.where(rhovals<0.9,(Ti0-Ti9)*(1-rhovals**2)**nu_Ti + Ti9,(Ti_sep-Ti9)/(0.1)*(rhovals-0.9)+Ti9)
ne     = np.where(rhovals<0.9,(ne0-ne9)*(1-rhovals**2)**nu_ne + ne9,(ne_sep-ne9)/(0.1)*(rhovals-0.9)+ne9)
Te     = np.where(rhovals<0.9,(Te0-Te9)*(1-rhovals**2)**nu_Te + Te9,(Te_sep-Te9)/(0.1)*(rhovals-0.9)+Te9)
nz     = rhovals**6 * ni * impfrac
#ni = np.linspace(niped,nisep,numpts)
#Ti = np.linspace(Tiped,Tisep,numpts)
#ne = np.linspace(neped,nesep,numpts)
#Te = np.linspace(Teped,Tesep,numpts)

####################################################################
####################################################################
## SPECIFY IMPURITY OR IMPURITIES
####################################################################
####################################################################

#        2:('Helium'),
#        3:('Lithium'),
#        4:('Berylium'),
#        5:('Boron'),
#        6:('Carbon'),
#        7:('Nitrogen'),
#        8:('Oxygen'),
#        10:('Neon'),
#        18:('Argon'),
#        36:('Krypton'),
#        54:('Xenon'),
#        86:('Radon')            


inucz = 86

p_plot = np.zeros(100)
p_plot_labels = []
#for inucz in [8,10,18,36,54, 86]:
for inucz in [6]:
    ####################################################################
    ####################################################################
    ## CREATE ADPAK INPUT FILE AND RUN ADPAK
    ####################################################################
    ####################################################################
    
    f = open('/home/max/Nextcloud/Max/max_phd/Codes/adpak/toadpak','w')
    
    f.write(' &inp')
    f.write('\n' + '  inucz = ' + str(inucz))
    f.write('\n' + '  zte = 2.0')
    f.write('\n' + '  zne = 2.0e14')
    f.write('\n' + '  laden = 1')
    f.write('\n' + '  ladtip = 1')
    f.write('\n' + '  leci = 1')
    f.write('\n' + '  ldrmlt = 2')
    f.write('\n' + '  ncxb = 0')
    f.write('\n' + '  ncxopt = 1')
    f.write('\n' + '  ivunit = 2')
    f.write('\n' + '  anneut = 1.0e11')
    f.write('\n' + '  vneut = 0.001')
    f.write('\n' + '  imode = 1')
    f.write('\n' + '  nte = 21')
    f.write('\n' + '  nne = 2')
    f.write('\n' + '  tei = 0.001 0.0015 0.003 0.005 0.01 0.02 0.03 0.05 0.07')
    f.write('\n' + '        0.1 0.2 0.3 0.4 0.5 0.7 1.0 2.0 3.0 5.0 10.0 50.0 100.0, 1000.0')
    f.write('\n' + '  anei = 1.0e+13 1.e+15')
    f.write('\n' + '  nmin = -3')
    f.write('\n' + '  nmax = 2')
    f.write('\n' + ' $end')
    f.write('\n' + '')
    f.write('\n')
    f.close()
    
    call(['./adpak'])
    #call(['/home/max/Nextcloud/Max/max_phd/Codes/ADPAK1/adpak'],cwd='/home/max/Nextcloud/Max/max_phd/Codes/ADPAK1/')
    
    ####################################################################
    ####################################################################
    ## READ IN ADPAK RESULTS
    ####################################################################
    ####################################################################
    
    
    #filepath = os.getcwd() + '/home/max/Nextcloud/Max_PhD/ADPAK1/outblk.dat'
    filepath = '/home/max/Nextcloud/Max/max_phd/Codes/adpak/outblk.dat'    
    with open(filepath, 'r') as line:
        dummy = next(line)
        data = re.findall(r'\d+', next(line))
        nte = int(data[0])
        nne = int(data[1])
        nchrgsr = int(data[2])
        
        ####################################################################
        ## Read in electron temperature section
        ####################################################################
        dummy = next(line)
        dummy = next(line)
    
        numlines = int(math.ceil(nte / 5))
        Te2 = []
        for i in range(0,numlines):
            data = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(line))
            for number in data:
                Te2.append(float(number))
        Tenp = np.asarray(Te2)
        T_max = np.amax(Tenp)
        T_min = np.amin(Tenp)
        Tenp = np.power(Tenp*0 + 10,Tenp)
    
        ####################################################################
        ## Read in electron density section
        ####################################################################
        dummy = next(line)
        dummy = next(line)
    
        numlines = int(math.ceil(nne / 5))
        ne2 = []
        for i in range(0,numlines):
            data = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(line))
            for number in data:
                ne2.append(float(number))
        
        ####################################################################
        ## Read in ionization rate section
        ####################################################################
        aliznr = []
        numlines = int(math.ceil(nte / 5))
        for j in range(0,nchrgsr):
            dummy = next(line)
            dummy = next(line)
        
            for i in range(0,numlines):
                data = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(line))
                for number in data:
                    aliznr.append(float(number))
        aliznrnp = np.reshape(np.asarray(aliznr),(-1,nte))
        aliznrnp = np.power(aliznrnp*0 + 10,aliznrnp)
        
        ####################################################################
        ## Read in radiation rate section
        ####################################################################
        alradr = []
        numlines = int(math.ceil(nne*nte / 5))
        for j in range(0,nchrgsr):
            dummy = next(line)
            dummy = next(line)
        
            for i in range(0,numlines):
                data = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(line))
                for number in data:
                    alradr.append(float(number))
        
        alradrnp = np.split(np.reshape(np.asarray(alradr),(-1,nne*nte)), 2, axis=1)[0] #(charge state,temperature,density) in that order
        alradrnp = np.power(alradrnp*0 + 10,alradrnp)
        if inucz==86:
            print (alradrnp[0].shape)
            print (alradrnp.shape)
            alradrnp = np.vstack((alradrnp[0],alradrnp))
    
        ####################################################################
        ## Read in recombination rate section
        ####################################################################
        alrecr = []
        numlines = int(math.ceil(nne*nte / 5))
        for j in range(0,nchrgsr):
            dummy = next(line)
            dummy = next(line)
        
            for i in range(0,numlines):
                data = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(line))
                for number in data:
                    alrecr.append(float(number))
        alrecrnp = np.split(np.reshape(np.asarray(alrecr),(-1,nne*nte)), 2, axis=1)[0] #(charge state,temperature,density) in that order
        alrecrnp = np.power(alrecrnp*0 + 10,alrecrnp)
        ####################################################################
        ####################################################################
        ## CORONAL EQUILBRIUM CALCULATION
        ####################################################################
        ####################################################################
        #Create nchrgsr x nchrgsr array
        M = np.zeros((nchrgsr,nchrgsr))
        
        #Create 2D interpolation functions
        inzr = interp2d(Tenp,np.arange(nchrgsr),aliznrnp,kind='linear')
        recr = interp2d(Tenp,np.arange(nchrgsr),alrecrnp,kind='linear') 
    
        Tpts = 100
        f = np.zeros((Tpts,nchrgsr))
        #f = np.zeros((numpts,nchrgsr))
        tcount = 0
        for T in np.logspace(T_min,1,Tpts):
        #for T in Ti:
            #print ("T = ",T)
            for index,value in np.ndenumerate(M):
                if index[0]==0: #do the first row differently
                    M[index[0],index[0]+0] = -(inzr(T,index[0]+0) + recr(T,index[0]+0)) 
                    M[index[0],index[0]+1] = recr(T,index[0]+1)
                elif index[0]==nchrgsr-1: #do the last row differently
                    M[index[0],index[0]-1] = inzr(T,index[0]-1)
                    M[index[0],index[0]+0] = -(inzr(T,index[0]+0) + recr(T,index[0]+0)) 
                else:
                    M[index[0],index[0]-1] = inzr(T,index[0]-1)
                    M[index[0],index[0]+0] = -(inzr(T,index[0]+0) + recr(T,index[0]+0)) 
                    M[index[0],index[0]+1] = recr(T,index[0]+1)
            M = M*1E10
            def null(M, eps=1e-20):
                u, s, vh = np.linalg.svd(M)
                #null_space = np.compress(min(abs(s)), vh, axis=0)
                spit_out = vh[np.argmin(np.abs(s))]
                #return null_space.T
                return spit_out
            
            f[tcount] = null(M).T #solve(M,b)
            tcount+=1
        f = np.where(f<0,-f,f)
        for row,vals in enumerate(f):
            f[row] = f[row]/np.sum(f[row])
            
        #make charge state fractional abundance interpolation
        f_int = interp2d(np.arange(nchrgsr),np.logspace(T_min,1,Tpts),f,kind='linear')
    ####################################################################    
    ####################################################################
    ## CALCULATE CHARGE STATE DENSITY PROFILES
    ####################################################################
    ####################################################################
    
    nz_matrix = np.zeros((numpts,nchrgsr))
    for i,rho in enumerate(rhovals):
        for j,cs in enumerate(np.arange(nchrgsr)):
            nz_matrix[i,j] = f_int(cs,Ti[i])*nz[i]
            
    
    ####################################################################    
    ####################################################################
    ## CALCULATE RADIATION PROFILES
    ####################################################################
    ####################################################################
    print (Tenp.shape, np.arange(nchrgsr).shape, alradrnp.shape)
    print ("Tenp = ",Tenp)
    print ("")
    print ("np.arange(nchrgsr) = ",np.arange(nchrgsr))
    print ("")
    print ("alradrnp = ",alradrnp)
    rad_matrix = np.zeros((numpts,nchrgsr))
    radr = interp2d(Tenp,np.arange(nchrgsr),alradrnp,kind='linear')
    for i,rho in enumerate(rhovals):
        for j,cs in enumerate(np.arange(nchrgsr)):
            rad_matrix[i,j] = radr(Ti[i],cs)
    
    rad_T_matrix = np.zeros((Tpts,nchrgsr))
    for i,tval in enumerate(np.logspace(T_min,1,Tpts)):
        for j,cs in enumerate(np.arange(nchrgsr)):
            rad_T_matrix[i,j] = radr(tval,cs)
            
    rrxf_T = np.sum(rad_T_matrix*f,axis=1)
            
    ####################################################################    
    ####################################################################
    ## CALCULATE VOLUME FUNCTION
    ####################################################################
    ####################################################################
    vol = np.zeros(rhovals.shape)
    for i,rho in enumerate(rhovals):
        rho2 = rho
        rho1 = rhovals[i-1]
        if i>0:
            vol[i] = 2 * pi**2 * R0_a * kappa * (rho2**2 - rho1**2)   
    ####################################################################    
    ####################################################################
    ## PLOT STUFF
    ####################################################################
    #################################################################### 
    ## GET ATOM NAME
    name = {
            2:('Helium'),
            3:('Lithium'),
            4:('Berylium'),
            5:('Boron'),
            6:('Carbon'),
            7:('Nitrogen'),
            8:('Oxygen'),
            10:('Neon'),
            18:('Argon'),
            36:('Krypton'),
            54:('Xenon'),
            86:('Radon'),
            74:('Tungsten')         
            } 
    
    abb = {
            2:('He'),
            3:('Li'),
            4:('Be'),
            5:('Bo'),
            6:('C'),
            7:('N'),
            8:('O'),
            10:('Ne'),
            18:('Ar'),
            36:('Kr'),
            54:('Xe'),
            86:('Ra'),
            74:('W')            
            } 
    
    plot_species = 0
    if plot_species == 1:   
        ## IONIZATION RATE BY CHARGE STATE
        fig1= plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.set_title(r'Ionization Rate by charge state',fontsize=20,y=1.08)
        #ax1.set_ylabel(r'$T_i\left(keV\right)$',fontsize=20)
        ax1.set_xlabel(r'Electron temperature $\left(keV\right)$',fontsize=15)
        ax1.grid(b=True,which='both',axis='both')
        ax1.set_xlim(10**T_min, 10**T_max)
        ax1.set_xscale('log',subsx=[-3,-2,-1,0,1,2,3])
        ax1.set_yscale('log',subsx=[-3,-2,-1,0,1,2,3])
        #ax1.set_ylim(0, 150)
        ax1.legend(loc='upper right', shadow=True)
        for i in range(0,nchrgsr):
            j = i/nchrgsr
            ax1.plot(Tenp,aliznrnp[i,:],label = i,color=plt.cm.jet(j))
        plt.tight_layout()
        
        ## RADIATION RATE BY CHARGE STATE
        fig2= plt.figure()
        ax1 = fig2.add_subplot(111)
        ax1.set_title(r'Radiation Rate by charge state', fontsize=20,y=1.08)
        #ax1.set_ylabel(r'$T_i\left(keV\right)$',fontsize=20)
        ax1.set_xlabel(r'Electron temperature $\left(keV\right)$',fontsize=15)
        ax1.grid(b=True,which='both',axis='both')
        ax1.set_xlim(10**T_min, 10**T_max)
        ax1.set_xscale('log',subsx=[-3,-2,-1,0,1,2,3])
        ax1.set_yscale('log',subsx=[-3,-2,-1,0,1,2,3])
        #ax1.set_ylim(0, 150)
        ax1.legend(loc='upper right', shadow=True)  
        for i in range(0,nchrgsr):
            j = i/nchrgsr
            plt.plot(Tenp,alradrnp[i,:],label = i,color=plt.cm.jet(j))
        plt.tight_layout()
        
        ## RECOMBINATION RATE BY CHARGE STATE
        fig3= plt.figure()
        ax1 = fig3.add_subplot(111)
        ax1.set_title(r'Recombination Rate by charge state',fontsize=20,y=1.08)
        #ax1.set_ylabel(r'$T_i\left(keV\right)$',fontsize=20)
        #ax1.set_xlabel(r'minor radius',fontsize=20)
        ax1.grid(b=True,which='both',axis='both')
        ax1.set_xlim(0, 1)
        #ax1.set_ylim(0, 150)
        ax1.legend(loc='upper right', shadow=True)        
        for i in range(0,nchrgsr):
            j = i/nchrgsr
            plt.plot(Tenp,alrecrnp[i,:],label=i,color=plt.cm.jet(j))
        plt.tight_layout()
        
        ## FRACTIONAL ABUNDANCES
        fig4= plt.figure()
        ax2 = fig4.add_subplot(111)
        ax2.set_title(r'%s fractional abundances'%(name[inucz]),fontsize=20,y=1.08)
        ax2.set_xlabel(r'Electron temperature $\left(keV\right)$',fontsize=15)
        ax2.set_ylabel(r'%s fractional abundances'%(name[inucz]),fontsize=15)
        ax2.grid(b=True,which='both',axis='both')
        ax2.set_xlim(10**T_min, 10**T_max)
        ax2.set_xscale('log',subsx=[-3,-2,-1,0,1,2,3])
        ax2.set_yscale('log',subsx=[-3,-2,-1,0,1,2,3])
        ax2.set_ylim(0.01, 1)
        for i in range(0,nchrgsr):
            j = i/nchrgsr
            ax2.plot(np.logspace(T_min,1,Tpts),f[:,i],label=r'$%s^{+%s}$'%(abb[inucz],i),color=plt.cm.jet(j),lw=3)
            #ax2.plot(Ti,f[:,i],label=r'$%s^{+%s}$'%(abb[inucz],i),color=plt.cm.jet(j),lw=3)
        for xc in [1,5]:
            ax2.axvline(x=xc,linewidth=3, color='black', ls='--')
        ax2.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=ceil(float(nchrgsr)/7), prop={'size':13.9}, shadow=True)  
        plt.tight_layout()
        
        ## CHARGE STATE DENSITY PROFILES
        fig5= plt.figure(figsize=(8, 6))
        ax1 = fig5.add_subplot(111)
        ax1.set_title(r'%s charge state density profiles'%(name[inucz]),fontsize=20,y=1.08)
        ax1.set_xlabel(r'normalized minor radius $\rho$',fontsize=15)
        ax1.set_ylabel(r'%s density $\left(m^{-3}\right)$'%(name[inucz]),fontsize=15)
        ax1.grid(b=True,which='both',axis='both')
        #ax1.set_xlim(10**T_min, 10)
        ax1.set_xlim(0, 1)
        ax1.set_ylim(0,0.6E17)
        #ax1.set_xscale('log',subsx=[-2,-1,0,1,2,3])
        #ax1.set_yscale('log',subsx=[-2,-1,0,1,2,3])
        #ax1.set_ylim(0.01, 1)
        for i in range(0,nchrgsr):
            j = i/nchrgsr
            ax1.plot(rhovals,nz_matrix[:,i],label=r'$%s^{+%s}$'%(abb[inucz],i),color=plt.cm.jet(j),lw=2)
        #ax1.plot(rhovals,ni,'--',color='black')
        #ax1.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=ceil(float(nchrgsr)/7), prop={'size':13.9}, shadow=True)  
        plt.tight_layout()
        
        ## RADIATION BY CHARGE STATE
        fig6= plt.figure(figsize=(8, 6))
        ax1 = fig6.add_subplot(111)
        ax1.set_title(r'%s radiation profiles by charge state'%(name[inucz]),fontsize=20,y=1.08)
        ax1.set_xlabel(r'normalized minor radius $\rho$',fontsize=15)
        ax1.set_ylabel(r'%s radiation $\left(W\right)$'%(name[inucz]),fontsize=15)
        ax1.grid(b=True,which='both',axis='both')
        #ax1.set_xlim(10**T_min, 10)
        ax1.set_xlim(0, 1)
        #ax1.set_ylim(0,0.6E17)
        #ax1.set_xscale('log',subsx=[-2,-1,0,1,2,3])
        #ax1.set_yscale('log',subsx=[-2,-1,0,1,2,3])
        #ax1.set_ylim(0.01, 1)
        for i in range(0,nchrgsr):
            j = i/nchrgsr
            ax1.plot(rhovals,nz_matrix[:,i]*rad_matrix[:,i]*ne*vol,label=r'$%s^{+%s}$'%(abb[inucz],i),color=plt.cm.jet(j),lw=2)
        #ax1.plot(rhovals,ni,'--',color='black')
        #ax1.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=ceil(float(nchrgsr)/7), prop={'size':13.9}, shadow=True)  
        plt.tight_layout()
        
        ## TOTAL RADIATION
        fig7= plt.figure(figsize=(8, 6))
        ax1 = fig7.add_subplot(111)
        ax1.set_title('%s total radiation profiles\nfor %s'%(name[inucz],machine),fontsize=20,y=1.08)
        ax1.set_xlabel(r'normalized minor radius $\rho$',fontsize=15)
        ax1.set_ylabel(r'%s radiation $\left(W\right)$'%(name[inucz]),fontsize=15)
        ax1.grid(b=True,which='both',axis='both')
        #ax1.set_xlim(10**T_min, 10)
        ax1.set_xlim(0, 1)
        #ax1.set_ylim(0,0.6E17)
        #ax1.set_xscale('log',subsx=[-2,-1,0,1,2,3])
        #ax1.set_yscale('log',subsx=[-2,-1,0,1,2,3])
        
        
        #ax1.set_ylim(0, 400)
        
        rad_tot_adpak = np.sum(np.sum(nz_matrix*rad_matrix,axis=1)*ne*vol)
        rad_tot_fit = np.sum(((1 + 0.3*Te) * 1E-43 * ne * nz * (nchrgsr-1)**(3.7-0.33*np.log(Ti)) * vol * 1E6)[:-1])
        ax1.text(0.05, 0.85, 
                 "Total ADPAK radiation = %s MW\nTotal fit radiation = %s MW"%(round(rad_tot_adpak/1E6,3),round(rad_tot_fit/1E6,3)), 
                 fontsize=15,
                 transform = ax1.transAxes,
                 bbox={'facecolor':'white', 'alpha':1, 'pad':10})
        ax1.text(0.05, 0.5, 
                 r'$n_z/n_i$ = %s '%(impfrac), 
                 fontsize=15,
                 transform = ax1.transAxes,
                 bbox={'facecolor':'white', 'alpha':1, 'pad':10})
        #plot adpak radiation profile
        ax1.plot(rhovals,np.sum(nz_matrix*rad_matrix,axis=1)*ne*vol,color='black',lw=2,label='ADPAK data')
        #plot 17.50 fit radiation profile
        ax1.plot(rhovals[:-1], ((1 + 0.3*Te) * 1E-43 * ne * nz * (nchrgsr-1)**(3.7-0.33*np.log(Ti)) * vol * 1E6)[:-1],'--',color='red',lw=2,label='Eq. 17.50')
        #print ("shapes = ",rhovals.shape,Te.shape)
        #ax1.plot(rhovals,Te,color='red',lw=2)
        #ax1.plot(rhovals,ni,'--',color='black')
        ax1.legend(bbox_to_anchor=(0.015, 0.6), 
                   loc='lower left', 
                   ncol=1, 
                   prop={'size':13.9},
                   #transform = ax1.transAxes,
                   shadow=True)  
        plt.tight_layout()
        
        ## DENSITY AND TEMPERATURE PROFILES
        fig8= plt.figure(figsize=(8, 6))
        ax1 = fig8.add_subplot(111)
        ax1.set_title('Electron Density and Ion Temperature Profiles\nfor %s'%(machine),fontsize=20,y=1.08)
        ax1.set_xlabel(r'normalized minor radius $\rho$',fontsize=15)
        ax1.set_ylabel(r'$n_e,n_z\left(m^{-3}\right)$',fontsize=15)
        ax1.grid(b=True,which='both',axis='both')
        #ax1.set_xlim(10**T_min, 10)
        ax1.set_xlim(0, 1)
        #ax1.set_ylim(0,0.6E17)
        #ax1.set_xscale('log',subsx=[-2,-1,0,1,2,3])
        #ax1.set_yscale('log',subsx=[-2,-1,0,1,2,3])
        #ax1.set_ylim(0.01, 1)
        
        ax1.plot(rhovals,ne,color='black',lw=2)
        ax1.plot(rhovals,nz,'--',color='black',lw=2)
        #ax1.plot(rhovals,ni,'--',color='black')
        #ax1.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=ceil(float(nchrgsr)/7), prop={'size':13.9}, shadow=True)  
        
        ax2 = ax1.twinx()
        #ax2.set_title(r'%s charge state radiation profiles'%(name[inucz]),fontsize=20,y=1.08)
        #ax2.set_xlabel(r'normalized minor radius $\rho$',fontsize=15)
        ax2.set_ylabel(r'Ion Temperature $\left(keV\right)$',fontsize=15,color='red')
        ax2.tick_params('y', colors='r')
        ax2.grid(b=True,which='both',axis='both')
        #ax2.set_xlim(10**T_min, 10)
        ax2.set_xlim(0, 1)
        #ax2.set_ylim(0,0.6E17)
        #ax2.set_xscale('log',subsx=[-2,-1,0,1,2,3])
        #ax2.set_yscale('log',subsx=[-2,-1,0,1,2,3])
        #ax2.set_ylim(0.01, 1)
        ax2.plot(rhovals,Te,'--',color='red',lw=2)
        plt.tight_layout()

    ## make array for plotting FIGURE 13.9 later
    p_plot = np.vstack((p_plot,rrxf_T))
    p_plot_labels.append(name[inucz])

p_plot = np.delete(p_plot, (0), axis=0)

plt_fig139 = 1
if plt_fig139 == 1:
    ## FIGURE 13.9
    fig9= plt.figure(figsize=(10, 6))
    ax1 = fig9.add_subplot(111)
    ax1.set_title(r'Impurity Radiative Power Loss',fontsize=25,y=1.08)
    ax1.set_xlabel(r'Electron temperature $\left(keV\right)$',fontsize=20)
    ax1.set_ylabel(r'$\frac{P}{n_e n_2}\left(W-m^3\right)$',fontsize=20)
    ax1.grid(b=True,which='both',axis='both')
    ax1.set_xlim(10**T_min, 10**T_max)
    ax1.set_xscale('log',subsx=[-3,-2,-1,0,1,2,3])
    ax1.set_yscale('log',subsx=[-3,-2,-1,0,1,2,3])
    #ax1.set_ylim(0.01, 1)
    for i,arr in enumerate(p_plot):
        ax1.plot(np.logspace(T_min,1,Tpts),arr,label=p_plot_labels[i], lw=2)
    ax1.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=1, prop={'size':13.9}, shadow=True)
    plt.tight_layout()