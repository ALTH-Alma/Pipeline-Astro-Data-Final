import os
import time 
import shutil
import argparse

def replace_sim(i):
    # Define the file name
    file_name = 'FARAD.py'

    # The new line to replace the existing line that starts with 'run'
    new_line = "run       = 'sim_{:d}' \n".format(i)    

    # Read the contents of the file
    with open(file_name, 'r') as file:
        lines = file.readlines()

    # Open the file in write mode to overwrite it
    with open(file_name, 'w') as file:
        for line in lines:
            # Check if the line starts with 'run'
            if line.startswith('run'):
                # Replace the line with the new line
                file.write(new_line)
            else:
                # Write the original line if it doesn't match
                file.write(line)

    print(f"The line starting with 'run' has been replaced with {new_line} in {file_name}.")


def replace_output(j):
    file_name = 'FARAD.py'
    new_line = "it        = {:d} \n".format(j)
    # Read the contents of the file
    with open(file_name, 'r') as file:
        lines = file.readlines()

    # Open the file in write mode to overwrite it
    with open(file_name, 'w') as file:
        for line in lines:
            # Check if the line starts with 'it'
            if line.startswith('it'):
                # Replace the line with the new line
                file.write(new_line)
            else:
                # Write the original line if it doesn't match
                file.write(line)

    print(f"The line starting with 'it' has been replaced with {new_line} in {file_name}.")


def replace_outdir(i,j):
    file_name = 'setup.py'
    new_line = "outdir = '/disk2/alma/pipeline/RT/sim_{:d}/output_{:d}' \n".format(i,j)
          
    with open(file_name, 'r') as file:
        lines = file.readlines()

    with open(file_name, 'w') as file:
        for line in lines:
            if line.startswith('outdir'):
                file.write(new_line)
            else:
                file.write(line)

    print(f"The outdir has been replced with {new_line} in {file_name}.")


def move_output_folder(folder, new_parent_folder):
    try:
        # Obtener la ruta completa del directorio que queremos mover
        full_path = os.path.abspath(folder)
        
        if os.path.isdir(full_path):
            # Obtener el nombre base del directorio
            folder_name = os.path.basename(full_path)
            # Construir la ruta completa de destino
            destination = os.path.join(new_parent_folder, folder_name)
            # Mover el directorio a la nueva ubicación
            shutil.move(full_path, destination)
            print(f"Se movió el directorio '{folder_name}' a '{destination}'.")
        else:
            print(f"Error: '{folder}' no existe o no es un directorio válido.")
    except (FileNotFoundError, OSError) as e:
        print(f"Error al intentar mover '{folder}': {e}")


int_time = time.time()  # Tiempo inicial en segundos desde época UNIX

N = 400  # Número de modelos de disco diferentes
n = 20    # Número de outputs por modelo
offset = 0  # Desfase en la numeración de los modelos
#new_folder = '/srv/nas02/alma/Radmc_data/out_leia_v2/All_sim_pix128/'
new_folder = '/srv/nas02/alma/Radmc_data/out_beam_v1/All_sim_pix128/'

for i in range(N):
    replace_sim(i + offset)
    for j in range(1, n):
        replace_output(j)
        replace_outdir(i + offset, j)
        os.system('python setup.py')
   
    output_folder_name = 'sim_' + str(i+offset)
    move_output_folder(output_folder_name, new_folder)
 

total_time = time.time() - int_time  # Tiempo total en segundos

# Imprime el tiempo total en formato legible
hours = int(total_time // 3600)
minutes = int((total_time % 3600) // 60)
seconds = int(total_time % 60)
print(f"Tiempo total de ejecución: {hours}h {minutes}m {seconds}s")
