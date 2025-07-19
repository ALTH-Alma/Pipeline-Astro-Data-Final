import subprocess
from multiprocessing import Pool

device = 1 
N = 96
offset = 304

def run_simulation(sim):
    new_par_file = 'setups/Alma_Vidal/Alma_Vidal_{:d}.par'.format(sim+offset)
    subprocess.run(["time", "./fargo3d", "-D", str(device), new_par_file])

if __name__ == '__main__':
    
    with Pool(40) as p:
       
        p.map(run_simulation, range(N))
