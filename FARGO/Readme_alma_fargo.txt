Para correr las simulaciones se deben seguir 4 pasos:

1- Generar configuraciones para las simulaciones con el script Alma_par_new.py modificando el número de simulaciones N y el offset, segun se necesite.
   (Esto genera archivos .par para las diferentes simulaciones que se quieran correr, asignando parametros de forma aleatoria, tiene configurado un time
    step de 0.31425 unidades de tiempo (2pi unidades de tiempo = 1 año orbital), Ntot = Ninterm, con Ntot aleatorio grande. Por lo que se tiene una salida
    con un tiempo de simulación largo)

2- Generar configuraciones de restart para las simulaciones con el script Alma_restart_par_new modificando el número de simulaciones N y el offset, segun se necesite.
   (Esto genera nuevos archivos de simulacion en base a las simulaciones anteriores, modificando solo el paso de tiempo a 0.0017799, Ntot a 200 e Ninterm 10, esto
    implica obtener cada 10 timesteps 1 salida hasta llegar a los 200 timestep (20 salidas)). Aproximadamente es como 1 salida por dias hasta llegar a los 20 dias.

3.- Correr las simulaciones 
