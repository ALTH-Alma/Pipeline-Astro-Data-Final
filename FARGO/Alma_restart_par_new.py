import os, fileinput, sys

N = 500 # Número de simulaciones
offset = 1500 # Offset de simulación

for sim in range(N):
    
    org_par_file = 'setups/Alma_Vidal_new/Alma_Vidal_new{:d}.par'.format(sim+offset)
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
