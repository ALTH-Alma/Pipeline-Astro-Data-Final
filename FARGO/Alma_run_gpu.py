### Corre simulaciones de forma secuencial, pero bajo el enfoque GPU de FARGO
### para make SETUP=Alma_Vidal_new PARALLEL=0 GPU=1.

### Para usarse antes se debe correr el 'make SETUP=Alma_Vidal_new PARALLEL=0 GPU=1'

import subprocess
import time

N = 100 # Número de simulaciones
offset = 0 # Offset de simulación

def run_simulation(sim):
    new_par_file = f'setups/Alma_Vidal_new/Alma_Vidal_new{sim}.par'
    subprocess.run(["time", "./fargo3d", new_par_file])

start = time.time()

for sim in range(N):
    sim_id = sim + offset
    print(f"Ejecutando simulación {sim_id}...")
    run_simulation(sim_id)

end = time.time()
print(f"Tiempo total: {end - start:.2f} segundos")
