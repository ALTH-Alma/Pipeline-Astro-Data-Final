Para correr las simulaciones FARGO3D se deben seguir 4 pasos:

1- Generar configuraciones para las simulaciones con el script 'Alma_par_new.py' modificando el número de simulaciones N y el offset, segun se necesite.
   (Esto genera archivos .par para las diferentes simulaciones que se quieran correr, asignando parametros de forma aleatoria, tiene configurado un time
    step de 0.31425 unidades de tiempo (2pi unidades de tiempo = 1 año orbital), Ntot = Ninterm, con un valor aleatorio grande de Ntot. Por lo que se tiene
    una salida con un tiempo de simulación largo).

Extra: Antes de correr las simulaciones se debe compilar el setup realizando un make, se puede hacer bajo 4 Enfoques. 
      - Secuencial CPU: make SETUP=Alma_Vidal_new PARALLEL=0 GPU=0
      - Secuencial GPU: make SETUP=Alma_Vidal_new PARALLEL=0 GPU=1
      - Paralelización interna (MPI) CPU: make SETUP=Alma_Vidal_new PARALLEL=1 GPU=0
      - Paralelización interna (MPI) GPU: make SETUP=Alma_Vidal_new PARALLEL=1 GPU=1

2.- Correr las simulaciones para las diferentes configuraciones .par utilizan los script 'Alma_run_X.py', donde X=[multiprocess, parallel, gpu] en base al enfoque
    seleccionado para correr las simulaciones, modificando el número de simulaciones (N), el offset y el número de nucleos para MPI (np), segun se necesite. 
    - Multiprocess es para 'secuencial CPU' pero generando multiples simulaciones en paralelo con Pool. 
    - Parallel es para usar con Paralelización interna (MPI) CPU.
    - Gpu es para usar con 'secuencial GPU'.

3- Generar configuraciones de restart para las simulaciones con el script 'Alma_restart_par_new.py' modificando el número de simulaciones N y el offset, segun se necesite.
   (Esto genera nuevos archivos de simulacion en base a las simulaciones anteriores, modificando en cada una, solo el paso de tiempo a 0.0017799, Ntot a 200 e Ninterm a 10,
    esto implica obtener cada 10 timesteps 1 salida hasta llegar a los 200 timestep (20 salidas en total)). Representando aproximadamente una salida por día durante 20 días
    de evolución simulada).

4.- Correr el restart de las simulaciones para las diferentes configuraciones .par que se modificaron utilizan los script 'Alma_run_restart_X.py', donde X=[multiprocess, 
    parallel, gpu] en base al enfoque seleccionado para correr las simulaciones, modificando el número de simulaciones (N), el offset y el número de nucleos para MPI (np),
    segun se necesite. La selección entre multiprocess, parallel y gpu funciona igual de con el paso 2. Si se desea cambiar el enfoque (por ejemplo, de CPU a GPU), se debe
    recompilar (make) con la opción correspondiente antes de ejecutar este paso.

Extra: La configuracion de la carpeta de outputs donde se quieran guardar las simulaciones se hace en el script 'Alma_par_new.py' al crear las configuraciones (paso 1).
