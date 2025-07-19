### Este script modifica los archivos de configuración de parámetros para reiniciar simulaciones.
### Copia archivos .par existentes y ajusta los valores de DT, Ninterm y Ntot
### para correr simulaciones más cortas a partir de un punto previo
### y obtener 20 salidas nuevas cada cierto intervalo de tiempo (según el valor de DT).

## DT: Es el tamaño del paso temporal (timestep) en las unidades de tiempo internas del código.
## 2π unidades de tiempo = 1 año orbital a 1 AU (unidad astronomica).

## Ntot: Índica cuánto tiempo dura la simulación en total (en pasos de tiempo).
## Nintem: Índica cada cuánto tiempo se guardan datos (en pasos de tiempo). 

### Cada 10 pasos de tiempo de ~0,0017799 unidades de tiempo (0.0103 días) se guarda una salida,
### hasta alcanzar los 200 pasos de tiempo (19,6 días). (10 pasos = 1.03 días)

import os, fileinput, sys

N = 500 # Número de simulaciones
offset = 1500 # Offset de simulación

for sim in range(N):
    
    org_par_file = 'setups/Alma_Vidal_new/Alma_Vidal_new_{:d}.par'.format(sim+offset)
    res_par_file = 'setups/Alma_Vidal_new/Alma_restart_new_{:d}.par'.format(sim+offset)
    os.system('cp '+org_par_file+' '+res_par_file)
    
    for line in fileinput.FileInput(res_par_file, inplace=True):
        if line.startswith("DT"):
            print("DT                      {:.6f}".format(0.00028328611*2*3.1415))
        elif line.startswith("Ninterm"):
            print("Ninterm                 {:d}".format(10))
        elif line.startswith("Ntot"):
            print("Ntot                    {:d}".format(200))
        else:
            sys.stdout.write(line)
    fileinput.close()
