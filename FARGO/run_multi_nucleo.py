import subprocess
import time

device = 1 
N = 2302 - 2139
offset = 2139

np = 30  # núcleos por simulación

# Lista de simulaciones ya completadas (parche por cambio de configuracion)
completadas = {
    2302, 2303, 2304, 2305, 2306, 2307, 2308, 2309, 2310, 2311, 2312, 2313,
    2500, 2501, 2502, 2698, 2896, 3094, 3095, 3292, 3293, 3294, 3295,
    3490, 3886, 3887, 3888
}

def run_simulation(sim):
    new_par_file = f'setups/Alma_Vidal_new/Alma_Vidal_new{sim + offset}.par'
    subprocess.run(["mpirun", "-np", str(np), "./fargo3d", new_par_file], check=True)

start = time.time()

for sim in range(N):
    sim_id = sim + offset
    if sim_id in completadas:
        print(f"Simulación {sim_id} ya completada, saltando...")
        continue

    print(f"Ejecutando simulación {sim_id}...")
    run_simulation(sim)

end = time.time()
print(f"Tiempo total: {end - start:.2f} segundos")
