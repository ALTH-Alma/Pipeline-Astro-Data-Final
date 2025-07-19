#### Script que crea archivos .par con diferentes configuraciones para simulaciones de discos
### usa definicion de parametros aleatoriamente con distribuciones
### crea configuraciones de planetas que se guardan en la carpeta planets

import fileinput, sys, random, shutil, subprocess, os
import numpy as np

# Parameters
N = 100 # Número de simulaciones
offset = 0 # Offset de la simulación
aspect_ratios = np.random.normal(0.05,0.025,N)
sigma0s = 10**np.random.normal(-3,0.5,N)
sigmaslopes = np.random.uniform(0.5,1.5,N)
flaringindex = np.random.uniform(0.0,0.5,N)
alphas = 10**np.random.uniform(-4,-2,N)
times = np.random.normal(0.5,0.25,N)*100.*1e-3/alphas
amax = np.random.uniform(0.1, 1, N)
rfac = np.random.uniform(2.0, 15.0, N)
ecc_big = np.random.uniform(0, 0.1, N)
ecc_small = np.random.uniform(0, 0.0001, N)
ymin = np.random.uniform(0.25, 0.5, N) 
nplanets = [1,1,1,1,2,2,3]
start_cfg = 1
planet_size = 1
org_par_file = 'setups/Alma_Vidal_new/Alma_Vidal_new.par'
device = 0

nps = []
mp1 = []
mp2 = []
mp3 = []
d1  = []
d2  = []
d3  = []
rp_out = []

# Writing planet file
for sim in range(N):
    Np = np.random.choice(nplanets)
    planetfile='planets/Alma_planets_new/{:d}.cfg'.format(sim+offset)
    dist1 = 1.0
    mass1 = 0.001*10**np.random.normal(0.0,1)
    first_planet = "Planet1         {:.2f}             {:.6f}    0.0            NO              NO\n".format(dist1,mass1)
    content = (
    "###########################################################\n"
    "#   Planetary system initial configuration\n"
    "###########################################################\n"
    "# Planet Name   Distance        Mass     Accretion      Feels Disk      Feels Others\n"
    f"{first_planet}"
    )
    
    # Write the content to the file
    with open(planetfile, 'w') as file:
        file.write(content)
    nps.append(Np)
    mp1.append(mass1)
    d1.append(dist1)
    if Np > 1:
        dist2 = np.random.uniform(1.5,2.5)
        mass2 = 0.001*10**np.random.normal(0.0,1)
        line_to_append = "Planet2         {:.2f}             {:.6f}    0.0            NO              NO".format(dist2,mass2)
        # Append the line to the file
        with open(planetfile, 'a') as file:
            file.write(line_to_append)
        mp2.append(mass2)
        d2.append(dist2)
    else:
        mp2.append(None)
        d2.append(None)
        rp_out.append(dist1)
    if Np > 2:
        dist3 = np.random.uniform(1.5*dist2,2*dist2)
        mass3 = 0.001*10**np.random.normal(0.0,1)
        line_to_append = "Planet3         {:.2f}             {:.6f}    0.0            NO              NO".format(dist3,mass3)
        # Append the line to the file
        with open(planetfile, 'a') as file:
            file.write(line_to_append)
        mp3.append(mass3)
        d3.append(dist3)
        rp_out.append(dist3)
    else:
        mp3.append(None)
        d3.append(None)
    if Np == 2:
        rp_out.append(dist2)
    print(f"The content has been written to {planetfile}")

print(nps)
print(mp1)
print(d1)
print(mp2)
print(d2)
print(mp3)
print(d3)
print(rp_out)


# Rewriting par file
for sim in range(N):
    new_par_file = 'setups/Alma_Vidal_new/Alma_Vidal_new_{:d}.par'.format(sim+offset)
    os.system('cp '+org_par_file+' '+new_par_file)
    ymax = np.random.uniform(2.0,3.0)*rp_out[sim]
    cutoff = np.random.uniform(1.5*rp_out[sim],1.5*ymax)
    ntot = int(times[sim]*20)

    for line in fileinput.FileInput(new_par_file, inplace=True):
        if line.startswith("AspectRatio"):
            print("AspectRatio             {:.2f}".format(aspect_ratios[sim]))
        elif line.startswith("Sigma0"):
            print("Sigma0                  {:.5f}".format(sigma0s[sim]))
        elif line.startswith("SigmaSlope"):
            print("SigmaSlope              {:.2f}".format(sigmaslopes[sim]))
        elif line.startswith("FlaringIndex"):
            print("FlaringIndex            {:.2f}".format(flaringindex[sim]))
        elif line.startswith("Alpha"):
            print("Alpha                   {:.5f}".format(alphas[sim]))
        elif line.startswith("rfactor"):
            print("rfactor                 {:.2f}".format(rfac[sim]))
        elif line.startswith("amax"):
            print("amax                    {:.2f}".format(amax[sim]))
        elif line.startswith("PlanetConfig"):
            print("PlanetConfig            planets/Alma_planets_new/{:d}.cfg".format(sim+offset))
        elif line.startswith("Eccentricity"):
            if planet_size == 0:
                print("Eccentricity            " + str(ecc_small[sim]))
            else:
                print("Eccentricity            " + str(ecc_big[sim]))
        elif line.startswith("Ymin"):
            print("Ymin                    {:.2f}".format(ymin[sim])) 
        elif line.startswith("Ymax"):
            print("Ymax                    {:.2f}".format(ymax))
        elif line.startswith("OutputDir"):
            print("OutputDir               @outputs/Alma_outputs_new/sim_{:d}".format(sim+offset))
        elif line.startswith("Cutoff"):
            print("Cutoff                  {:.2f}".format(cutoff))
        elif line.startswith("Ntot"):
            print("Ntot                    {:d}".format(ntot))
        elif line.startswith("Ninterm"):
            print("Ninterm                 {:d}".format(ntot))
        else:
            sys.stdout.write(line)
    fileinput.close()
