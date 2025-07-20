### Corre simulaciones de forma secuencial, pero bajo el enfoque parallelo de FARGO
### para make SETUP=Alma_Vidal_new PARALLEL=1 GPU=0
### Implica paralelizar internamente usando múltiples núcleos (MPI).

### Para usarse antes se debe correr el 'make SETUP=Alma_Vidal_new PARALLEL=1 GPU=0'

import subprocess
import time

N = 100 # Número de simulaciones
offset = 0 # Offset de simulación

np = 30  # Núcleos por simulación

def run_simulation(sim):
    intStart = time.time() 
    new_par_file = f'setups/Alma_Vidal_new/Alma_Vidal_new_{sim}.par'
    subprocess.run(["mpirun", "-np", str(np), "./fargo3d", new_par_file], check=True)
    intEnd = time.time() 
    rint(f"Simulación {sim+offset} completada en {intEnd - intStart:.2f} segundos")


start = time.time()

for sim in range(N):
    sim_id = sim + offset
    print(f"Ejecutando simulación {sim_id}...")
    run_simulation(sim_id)

end = time.time()
print(f"Tiempo total: {end - start:.2f} segundos")
