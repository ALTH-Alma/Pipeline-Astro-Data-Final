# coding: utf8
# This file summarizes some important parameters and setups.
import numpy as np
import sys
import os
nthreads = 20      # Number of threads used for MCMC
it        = 10 
dust      = 1     # Dust deduced from gas (0) or Fluid output for dust (1)
outstyle  = 1     # Output format of FARGO3D simulation. (1) for .dat files, (2) for .mpio files
N         = 4     # Number by which each dimension is reduced (fargo->radmc)
                  # to make it faster for radmc3d. I.e. radmc grid will be (nx//N,ny//N) 
dim       = 2     # Number of dimensions of FARGO3D run. 2 assumes an r,phi - grid
run       = 'sim_255' 
#indir     = '/disk2/alma/pipeline/hydro/outputs/Alma_outputs_new/'+run+'/'
#indir     = '/srv/nas02/alma/Fargo_data/outputs/out_leia_v2/'+run+'/'
indir     = '/srv/nas02/alma/Fargo_data/outputs/out_beam_v1/'+run+'/'

Sizes     = np.loadtxt(indir+"dust.par")[:,1] # dust sizes in cm multidust config 1
#Sizes     = np.loadtxt(indir+"dust.par")[1:] # dust sizes in cm unidust config 1
 
if it == 1:
    os.system('mkdir '+run)
os.system('mkdir '+run+'/output_{:d}'.format(it))


file_name = indir+'variables.par'

# Initialize the variable to store the value of rfactor
rfactor_value = None

# Open the file and read its contents
with open(file_name, 'r') as file:
    for line in file:
        # Check if the line starts with 'rfactor'
        if line.startswith('RFACTOR'):
            # Extract the value after 'rfactor' and any whitespace
            rfactor_value = float(line.split()[1])
            break  # Exit the loop once the value is found

# Print the extracted value
print(f"The value of rfactor is: {rfactor_value}")
R_Scaling = rfactor_value*5.2



micro     = [0.0001] # This is to add very small dust 

add_micro = False

if(add_micro==True):
    Sizes = np.concatenate([micro,Sizes])
ndust = len(Sizes) # Number of dust fluids                          
dsi = 0          # dust spacing interpolation scheme
innerdisk = False # Set this to True if there should be an artifical inner disc implemented (can be ignored for now)
if(innerdisk == False):
    N_innerdisk = 0
    N_outerdisk = 0
    CutEdges = True  # if True, the dampingzones of the simulation (which could hold unwanted structure) are cut off 
else:
    CutEdges = False
    N_innerdisk = 128
    N_outerdisk = 0
noemission = 0.05   # Radius, unto where the emission is not considered in the ray-tracing
quality = "low"     # set quality to 'low' or 'high'. This set number of cells and number of photons. 
if (quality == "high"):
    N=2
rotate_phi=float(-0.0)  # We can rotate the simulation by some angle (here, in degrees)
