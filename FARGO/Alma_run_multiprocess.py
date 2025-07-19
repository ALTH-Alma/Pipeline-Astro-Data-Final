### Deja corriendo multiples simulaciones de fargo en paralelo pero en el enfoque secuencial de fargo
### para make SETUP=Alma_Vidal_new PARALLEL=0 GPU=0

### Antes de usarse se debe correr el 'make SETUP=Alma_Vidal_new PARALLEL=0 GPU=0'

import subprocess
from multiprocessing import Pool

N = 96 # Número de simulaciones
offset = 0 # Offset de simulación

def run_simulation(sim):
    new_par_file = 'setups/Alma_Vidal/Alma_Vidal_{:d}.par'.format(sim+offset)
    subprocess.run(["time", "./fargo3d", new_par_file])

if __name__ == '__main__':
    
    with Pool(40) as p:
        p.map(run_simulation, range(N))
