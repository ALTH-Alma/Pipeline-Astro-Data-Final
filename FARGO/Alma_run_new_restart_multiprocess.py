import subprocess
from multiprocessing import Pool

device = 0
N = 500
offset = 1500
restart_file = 1

pool = 100

def run_simulation(sim):
    new_par_file = 'setups/Alma_Vidal_new/Alma_restart_new_{:d}.par'.format(sim+offset)
    subprocess.run(["time", "./fargo3d", "-D", str(device), "-S", str(restart_file), new_par_file])

if __name__ == '__main__':
    
    with Pool(pool) as p:
       
        p.map(run_simulation, range(N))
