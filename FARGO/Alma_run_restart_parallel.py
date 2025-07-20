## Corre el restart de la configuración en paralelo de FARGO

import subprocess
import time

N = 100 # Número de simulaciones
offset = 0 # Offset de simulación
restart_file = 1 # Archivo de salida de la simulacion desde el cual se inicia el restart

np = 30  # Núcleos por simulación

def run_simulation(sim):
    new_par_file = 'setups/Alma_Vidal_new/Alma_restart_new_{:d}.par'.format(sim)
    subprocess.run(["mpirun", "-np", str(np), "./fargo3d", "-S", str(restart_file), new_par_file], check=True)

start = time.time()

for sim in range(N):
    sim_id = sim + offset
    print(f"Ejecutando simulación {sim_id}...")
    run_simulation(sim_id)

end = time.time()
print(f"Tiempo total: {end - start:.2f} segundos")
