import subprocess

device = 0 
N = 13
offset = 287
restart_file = 1
for sim in range(N):
    new_par_file = 'setups/Alma_Vidal/Alma_restart_{:d}.par'.format(sim+offset)
    #Run simulation
    subprocess.run(["time", "./fargo3d", "-D", str(device), "-S", str(restart_file),new_par_file])
