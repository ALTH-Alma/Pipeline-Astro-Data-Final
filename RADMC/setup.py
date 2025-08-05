# This file holds the main flow of the RT pipeline.
import numpy as np
import math 
import copy
import sys
from functions import *
from Gausssmooth import *
from FARAD import *
from interinner import *

outdir = '/disk2/alma/pipeline/RT/sim_255/output_16' 

### revisar flerring shadows
# -------------------------------------------------------------------  
# set disk Model
# -------------------------------------------------------------------  
print("Defining Disk model")
M = Model()                       # radmc 

# -------------------------------------------------------------------  
# write setup files
# -------------------------------------------------------------------  
mesh   = Mesh(directory=indir)
params = Parameters(directory=indir)
sigmag = Field_dat(mesh,field='gas', quantity='dens',directory=indir)
'''
In the isothermal version of Fargo3d, the output quantity named 'energy'
contains the sounds speed cs. To calculate the temperature from that,
one can use the formula T = cs^2*µ*mp/k_b, which will be done when writing
the file.
'''
T      = Field_dat(mesh,field='gas', quantity='energy',directory=indir)
radtrans     = True
Writedust    = True
WriteOpac    = True
Dototalint   = True
DoRayTrace   = True

print("writing AMR grid")
write_AMRgrid(mesh,params)    
print("writing wavelength")
wavegrid=write_wavelength(M)

if WriteOpac:
    for i in range(1,ndust+1):

        os.system('optool pyr -s -radmc -l wavelength_micron.inp -amin {:.5f} -amax {:.5f} -na 50'.format(Sizes[i-1]*1.0e4/np.sqrt(10.0), Sizes[i-1]*1.0e4*np.sqrt(10.0)))
        os.system('mv dustkapscatmat.inp dustkapscatmat_{:d}.inp'.format(i))
    
    with open('dustopac.inp','w+') as f:
        f.write('2               Format number of this file\n')
        f.write('{:d}               Nr of dust species\n'.format(ndust))
        f.write('============================================================================\n')
        for i in range(1,ndust+1):
            f.write('10               Way in which this dust species is read\n')
            f.write('0               0=Thermal grain\n')
            f.write('{:d}        Extension of name of dustkappa_***.inp file\n'.format(i))
            f.write('----------------------------------------------------------------------------\n')

if Writedust:
    '''
    The dust species from Fargo3d a discrete - they represen one single size.
    To get a more continuous distribution for the opacities, we let the fargo density represent
    a range of dust sizes, rather than a single size. We then calculate opacities using 'optool'.
    Source files and documentation can be found at https://github.com/cdominik/optool.
    '''
    
    '''
    We load in the dust densities from FARGO3D. These are surface densities in the cas of 2D runs.
    They are in code units. Therefore, we need to convert them into 3d densities by using an
    analytical function and we have to multiply them by their cgs units. In the case of adding
    µm dust grains, we can take the gas surface denisty multiplied by some factor epsilon.
    '''
    sigmad = []
    if(add_micro==True):
        for i in range(1,ndust+1):
            if(i<len(micro)+1):
                sigmad.append(Field_dat(mesh,field='gas',quantity='dens',directory=indir))
                sigmad[i-1] *= 0.0001
            
            else:
                sigmad.append(Field_dat(mesh,field='dust'+str(i-len(micro)), quantity="dens",directory=indir))
            write_dust_densities(sigmad[i-1],sigmag, M, params, mesh, i, Sizes[i-1],ndust)
            write_dust_temp(T,params, mesh,i,ndust)
    else:
        for i in range(1,ndust+1):
            sigmad.append(Field_dat(mesh,field='dust'+str(i), quantity="dens",directory=indir))
            write_dust_densities(sigmad[i-1], sigmag, M, params, mesh, i, Sizes[i-1],ndust)
            write_dust_temp(T,params, mesh,i,ndust)
       
print("writing radmc3d input")
if quality == 'high':
    write_radmc3dinp(nphot=1e8, nphot_scat=1e7, nphot_spec = 1e6, modified_random_walk=1, scattering_mode_max=2, setthreads=nthreads)
else:
    write_radmc3dinp(nphot=1e7, nphot_scat=1e6, nphot_spec = 1e5, modified_random_walk=1, scattering_mode_max=2, setthreads=nthreads)
run_raytracing(M)
exportfits(M, Plot=False)
os.system('python run_convol_ALMA.py')
os.system('mv *.fits *.out '+outdir+'/.')
os.system('rm *.inp')
