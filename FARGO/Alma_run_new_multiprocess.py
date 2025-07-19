import subprocess
from multiprocessing import Pool
import time

device = 1 
N = 7900
offset = 2104

parallel = True  # Cambiar esto según el tipo de paralelismo que quieras
np = 10
pool = 10

#Parallel setup
#make SETUP=Alma_Vidal_new PARALLEL=1 GPU=0

#Secuencial GPU setup
#make SETUP=Alma_Vidal_new PARALLEL=0 GPU=1

def run_simulation(sim):
    new_par_file = 'setups/Alma_Vidal_new/Alma_Vidal_new{:d}.par'.format(sim + offset)

    # Configuración Parallel CPU
    if parallel:
        # Usa `mpirun` para ejecutar en paralelo
        subprocess.run(["time", "mpirun", "-np", str(np), "./fargo3d", new_par_file])
    else:
        # Configuración Sequential GPU
        subprocess.run(["time", "./fargo3d", new_par_file])


# Actualmente corriendo con parallel cpu
if __name__ == '__main__':
    start = time.time()

    with Pool(pool) as p:
        p.map(run_simulation, range(N))

    end = time.time()

    print(f"Tiempo total: {end - start:.2f} segundos")

