# coding: utf8
import numpy as np
import sys
import os
from FARAD import *

def extend_inner_disk_rad():
    r_in = []
    r_out = []
    for i in range(N_innerdisk):
        if system == "V4046" and Artificial_Planet==True:
            in_rad = 0.7*np.sqrt(np.loadtxt(indir+'planet_art.dat')[1]**2+np.loadtxt(indir+'planet_art.dat')[2]**2)
            out_rad = min(1.3*np.sqrt(np.loadtxt(indir+'planet_art.dat')[1]**2+np.loadtxt(indir+'planet_art.dat')[2]**2),0.4)
        else:
            in_rad = 0.01
            out_rad= 0.1
        r_in.append(in_rad*(out_rad/in_rad)**(float(i)/(1.0*N_innerdisk)))
            
    for i in range(N_outerdisk):
        r_out.append(4.0*(6.5/4.0)**(float(i)/(1.0*N_outerdisk)))
    r_old = np.loadtxt(indir+"domain_y.dat")[3:-3]
    r_update = np.concatenate([r_in,r_old])
    r_new = np.concatenate([r_update,r_out])
    np.savetxt('radial_domain.txt', r_new)
    return r_new

def extend_inner_disk(r,dens_old,fluid_name):
    au = 14959787066000.0           # [Astronomical unit] = cm
    dens_in = np.array(0.0*dens_old[:N_innerdisk,:])
    dens_out= np.array(0.0*dens_old[:N_outerdisk,:])
    dens_step = np.concatenate([dens_in,dens_old],axis=0)
    dens = np.concatenate([dens_step,dens_out],axis=0)
    Sigma_c = 1000.
    Sigma_c *= (R_Scaling*au)**2/(1.75*2e33)
    R_c = 16./R_Scaling
    Rsub = 0.03#/R_Scaling
    #Rin = 5./R_Scaling
    Rin = 0.5#3./R_Scaling
    Rincav = 9.0/R_Scaling
    Ractive = 8.0/R_Scaling
    delta_gas_in = 5.0e-6
    delta_gas_cavinner = 5.0e-3
    gamma = 0.5
    #ep1 = 0.0005
    ep1 = 0.00001
    '''
    print R_Scaling
    if add_micro == True:
        name = 'micro'
    else:
        name = 'dust1'
    newradius = open("newrad.txt",'w')
    for ir in range(np.shape(dens)[0]):
        newradius.write(str(ir)+'\t'+str(r[ir]) + '\n')
        for ip in range(np.shape(dens)[1]):
            if (r[ir] < 0.0003):
                dens[ir,ip] = 0.0     
            elif (Rsub <= r[ir] < Rin):  # inner part
                #if fluid_name == 'gas' or fluid_name == name:
                dens[ir,ip] = Sigma_c  * (r[ir]/0.5)**(-gamma)
                dens[ir,ip] *= np.exp(-(r[ir]/R_c)**(2.0-gamma))
                dens[ir,ip] *= delta_gas_in
                if fluid_name == name:
                    dens[ir,ip] *= 0.0*ep1
                if fluid_name == 'dust2':
                    dens[ir,ip] *= 0.0*0.01
                
            elif (Rin <= r[ir] < Rincav):  # first disk
                #if fluid_name == name or fluid_name == 'gas':
                dens[ir,ip] = Sigma_c  * (r[ir]/R_c)**(-gamma)
                dens[ir,ip] *= np.exp(-(r[ir]/R_c)**(2.-gamma))
                dens[ir,ip] *= delta_gas_cavinner
                if fluid_name == name:
                    dens[ir,ip] *= ep1
                if fluid_name == 'dust2' or fluid_name == 'dust3':
                    dens[ir,ip] *= 0.0
            elif (r[ir]>=Rincav):
                if fluid_name == 'micro':
                    dens[ir,ip] = 0.0
            #if r[ir]<=Ractive: and fluid_name != 'gas' and fluid_name != name:
            #    dens[ir, ip] = 0.0
            if r[ir] >= 4.5:
                dens[ir,ip] = 0.0
    newradius.close()
    '''
    return dens


