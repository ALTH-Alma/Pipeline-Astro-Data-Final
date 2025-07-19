### Deja corriendo multiples simulaciones de fargo en paralelo pero en el enfoque secuencial de fargo
### para make SETUP=fargo PARALLEL=0 GPU=0

import subprocess
from multiprocessing import Pool

N = 96
offset = 304

def run_simulation(sim):
    new_par_file = 'setups/Alma_Vidal/Alma_Vidal_{:d}.par'.format(sim+offset)
    subprocess.run(["time", "./fargo3d", new_par_file])

if __name__ == '__main__':
    
    with Pool(40) as p:
       
        p.map(run_simulation, range(N))
