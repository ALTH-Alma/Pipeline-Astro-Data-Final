import subprocess
import time

device = 1 
N = 7900
offset = 2107

np = 8  # núcleos por simulación
max_concurrent = 12  # máximo de simulaciones corriendo al mismo tiempo

# Lista de simulaciones ya completadas (parche por cambio de configuracion)
completadas = {
    2302, 2303, 2304, 2305, 2306, 2307, 2308, 2309, 2310, 2311, 2312, 2313,
    2500, 2501, 2502, 2698, 2896, 3094, 3095, 3292, 3293, 3294, 3295,
    3490, 3886, 3887, 3888
}

def run_simulation(sim):
    new_par_file = f'setups/Alma_Vidal_new/Alma_Vidal_new{sim + offset}.par'
    proc = subprocess.Popen(["mpirun", "-np", str(np), "./fargo3d", new_par_file])
    return proc

start = time.time()
running_procs = []

for sim in range(N):
    sim_id = sim + offset
    if sim_id in completadas:
        print(f"Simulación {sim_id} ya completada, saltando...")
        continue

    proc = run_simulation(sim)
    running_procs.append(proc)

    while len(running_procs) >= max_concurrent:
        for p in running_procs:
            if p.poll() is not None:
                running_procs.remove(p)
                break
        else:
            time.sleep(1)

# Esperar a que terminen todos
for p in running_procs:
    p.wait()

end = time.time()
print(f"Tiempo total: {end - start:.2f} segundos")
