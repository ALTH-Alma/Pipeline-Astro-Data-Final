# Corre el restart de las simulaciones para la configuración miltiprocess

import subprocess

N = 100 # Número de simulaciones
offset = 0 # Offset de simulación
restart_file = 1 # Archivo salida de la simulacion desde la que se empieza el restart (desde la primera)

for sim in range(N):
    new_par_file = 'setups/Alma_Vidal_new/Alma_restart_new_{:d}.par'.format(sim+offset)
    subprocess.run(["time", "./fargo3d", "-S", str(restart_file), new_par_file])
