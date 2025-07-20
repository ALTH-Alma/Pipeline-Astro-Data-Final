# Corre el restart de las simulaciones para la configuración miltiprocess

import subprocess
from multiprocessing import Pool

N = 100 # Número de simulaciones
offset = 0 # Offset de simulación
restart_file = 1 # Archivo salida de la simulacion desde la que se empieza el restart (desde la primera)

def run_simulation(sim):
    new_par_file = 'setups/Alma_Vidal_new/Alma_restart_new_{:d}.par'.format(sim+offset)
    subprocess.run(["time", "./fargo3d", "-S", str(restart_file), new_par_file])

if __name__ == '__main__':
    
    with Pool(40) as p:
        p.map(run_simulation, range(N))
