import numpy as np
import math as ma
import copy
import astropy.io.fits as pyfits
import os
import sys
import re
from FARAD import *
from interinner import *
from math import erf
from math import atan2

# =================================================================== 
#                       Parametric Disk Model
#            by Seba Perez, Simon Casassus and Seba Marino
# =================================================================== 
# 
# Considerations: 
#    - Dust is assumed standard dust size distribution.
# 
#    - Spherical coordinate system RADMC3D
# 
#        x = r     == radius
#        y = theta == colatitude 
#        z = phi   == azimuthal angle 
#

# For input and output data to file, for stuff on a regular grid, the
# order of nested loops over coordinates would be:

#       do iz=1,amr_grid_nz
#          do iy=1,amr_grid_ny
#             do ix=1,amr_grid_nx
#                << read or write your data >>
#             enddo
#          enddo
#       enddo
# 
# in RADMC3D
# x, y, z <---> r, theta, phi
# 
# 
# 
# =================================================================== 



# -------------------------------------------------------------------  
# global constants CGS
# -------------------------------------------------------------------  
M_Sun = 1.9891e33               # [M_sun] = g
G = 6.67259e-8                  # [G] = dyn cm2 g-2 
au = 14959787066000.0           # [Astronomical unit] = cm
m_H = 1.673534e-24              # [Hydrogen mass] = g
pi = 3.141592653589793116       # [pi]
R_Sun = 6.961e10                # [R_sun] = cm
kB = 1.380658E-16
pc = 3.085678E18                # [pc] = cm
c  = 2.99792458E10              # [c] = cm s-1
sigmaB=5.6704E10-5  # erg cm-2 s-1 K-4.

# -------------------------------------------------------------------
# Class "Model" defines a parametric disk model
# -------------------------------------------------------------------

# Aca se cambia la incl (inclinacion) y la posang (rotacion)
class Model():
    def __init__(self, nrad=128, nsec=256, ncol=16, # defaults
                 wmin=0.09,
                 wmax=3000.0,
                 Nw=100,
                 distance = 100,     # distance in parsec
                 Lambda=1300, # observing wavelength in micrometre
                 #npix=1024, # number of pixels for the image
                 #npix=300,
                 npix=128,
                 incl=0.0,
                 posang=0.0):
        # grid pars
        initlocals=locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            #if VerboseInit:
            #    print( "DConeMaps setting ",a_attribute," to ",initlocals[a_attribute])
            setattr(self,a_attribute,initlocals[a_attribute])

        self.nrad = nrad
        self.nsec = nsec
        self.ncol = ncol
        self.nx = nrad
        self.ny = ncol
        self.nz = nsec
        self.distance = distance * pc
        self.Lambda = Lambda
        self.npix   = npix
        #self.npix = cal_pix(M, R_Scaling)
        self.incl   = incl
        self.posang   = posang
        self.nx = nrad
        self.ny = ncol
        self.nz = nsec

# -------------------------------------------------------------------
# building mesh arrays for r, theta, phi (x, y, z)
# -------------------------------------------------------------------  
# write AMR Grid
# -------------------------------------------------------------------  
# -------------------------------------------------------------------  
# building mesh arrays for theta, r, phi (x, y, z)
# -------------------------------------------------------------------  



class Mesh():
    # based on P. Benitex routine
    """
    Mesh class, keeps all mesh data.
    Input: directory [string] -> place where domain files are
    """
    def __init__(self,directory="",staggered='c'):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        try:
            domain_x = np.loadtxt(directory+"domain_x.dat") # AZIMUTH
        except IOError:
            print("IOError with domain_x.dat")
        try:
            # avoid ghost cells!
            domain_y = np.loadtxt(directory+"domain_y.dat")[3:-3] # RADIUS
            if(innerdisk == True):
                domain_y = extend_inner_disk_rad()
        except IOError:
            print("IOError with domain_y.dat")
        try:
            # avoid ghost cells!
            domain_z = np.loadtxt(directory+"domain_z.dat")[3:-3] # COLATITUDE
        except:
            print("No vertical domain found - Is the simulation 2D?")

        # --------------------------------------------------------------
        # SWAPPING CONVENTION
        # --------------------------------------------------------------
        # rad will be stored in X
        # col will be stored in Y
        # azi will be stored in Z
        # X, Y, Z --> Y, Z, X
        # --------------------------------------------------------------
        self.x_full = 0.5*(domain_y[:-1] + domain_y[1:]) # Needed to compute density field in extend_inner()
        domain_x = domain_x[::N]
        domain_y = domain_y[::N]
        self.xm = domain_y # rad-edge
        if (dim == 3):
            self.ym = domain_z # col-Edge
        else:
            thmin = np.pi/2. 
            thmax = np.pi/2. + np.pi/4.0
            
            dth_min = 0.004
            ny = 32
            xi = np.log((thmax-thmin)/dth_min)/np.log(17.0)
            ym = []
            for j in range(ny//2):
                ym.append(thmin + dth_min*j**(xi))
            ym = np.array(ym)
            ymud = np.pi - ym
            ymud = ymud[::-1]
            ym = np.concatenate([ymud[:-1],ym])

            #thmin = pi/2. - pi/8.
            #thmax = pi/2.
            #ymp = np.linspace(np.log10(thmin),np.log10(thmax),17)
            #ym = -1.0*10**(ymp)+thmin+thmax # finer towards the midplane
            #ymud = pi - ym
            #ym = ym[::-1]
            #ym = np.concatenate([ym,ymud[1:]])
            self.ym = ym
        self.zm = domain_x # azi-Edge

        self.xmed = 0.5*(domain_y[:-1] + domain_y[1:]) # rad-Center
        if (dim == 3):
            self.ymed = 0.5*(domain_z[:-1] + domain_z[1:]) # col-Center
        else:
            self.ymed = 0.5*(ym[:-1] + ym[1:]) # Y-Center theta                                                                                             
        self.zmed = 0.5*(domain_x[:-1] + domain_x[1:]) # azi-Center
        
        # surfaces taken from the edges
        # make 2D arrays for x, y, that are (theta, r)
        T,R = np.meshgrid(self.zm, self.xm)
        R2  = R*R

        # --------------------------------------------------------------
        # plotting the mesh
        # --------------------------------------------------------------

        # surfaces taken from the edges
        # make 2D arrays for rad, azi (r, phi)
        # radius vs azimuth (cartesian)
        Plot = False
        if Plot:
            import matplotlib.pyplot as plt
            ax = plt.gca()
            P,R = np.meshgrid(self.zm, self.xm)    
            X = R*np.cos(P)
            Y = R*np.sin(P) 
            plt.pcolor(X,Y,np.random.rand(len(self.xm),len(self.zm)), 
                       cmap='plasma', edgecolors='black')
            plt.axis('equal')
            plt.show()

        self.surf = 0.5*(T[:-1,1:]-T[:-1,:-1])*(R2[1:,:-1]-R2[:-1,:-1])
        # fixing grid sizes due to reordering
        self.nx = (self.xm.size-1)
        #if(innerdisk==True):
        #    self.nx += 200
        self.ny = (self.ym.size-1)
        self.nz = (self.zm.size-1)

        # --------------------------------------------------------------
        # Unit conversion for radii 
        # assumes R0 = 1 au
        # --------------------------------------------------------------
        self.xm *= R_Scaling*au
        self.xmed *= R_Scaling*au

        # now, staggering:
        if staggered.count('x')>0:
            self.x = self.xm[:-1] # do not dump last element
        else:
            self.x = self.xmed
        if staggered.count('y')>0:
            self.y = self.ym[:-1]
        else:
            self.y = self.ymed
        

# -------------------------------------------------------------------  
# reading parameter file
# -------------------------------------------------------------------  

class Parameters():
    # based on P. Benitex routine
    """
    Reading simulation parameters.
    input: string -> name of the parfile, normally variables.par
    """
    def __init__(self, directory=''):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        try:
            params = open(directory+"variables.par",'r') # opening parfile
        except IOError:         # error checker.
            print(paramfile + " not found.")
            return
        lines = params.readlines()     # reading parfile
        params.close()                 # closing parfile
        par = {}                       # allocating a dictionary
        for line in lines:             # iterating over parfile
            name, value = line.split() # spliting name and value (first blank)
            try:
                float(value)           # first trying with float
            except ValueError:         # if it is not float
                try:
                    int(value)         # we try with integer
                except ValueError:     # if it is not integer, we know it is string
                    value = '"' + value + '"'
            par[name] = value          # filling variable
        self._params = par             # control atribute, good for debbuging
        for name in par:               # iterating over the dictionary
            exec("self."+name.lower()+"="+par[name]) # making atributes at runtime
        # Now we read in the dust sizes:
        try: 
            dsizes = open(directory+"dust.par")
        except IOError:
            print("dust.par not found!")
            return
        dlines = dsizes.readlines()
        dsizes.close()
        size = {}
        for dline in dlines:
            number, dsize = line.split()
            size[number] = dsize
        for number in size:
            exec("self.size"+str(number)+"="+size[number])
        #This saves the dust sizes as self.size1, self.size2, etc

# Get the pixel
def cal_pix(M, R_Scaling):
    G,au = 6.67e-8, 1.5e13
    R_cm = R_Scaling*au
    vK = np.sqrt(G*M.Mstar/R_cm)
    #How far does planet move in one year?
    yr = 3.154e7
    ds = vK*yr
    #How many arcsec does that correspond to?
    ds_au = ds/au
    ds_arcsec = ds_au/M.dist
    LAS = 3.0
    factor = 3
    npix = LAS/ds_arcsec*factor
    return npix

#### Calling and reshaping the fields
def Field_dat( p, field, quantity='',staggered='c', directory='', dtype='float64'):
    if len(directory) > 1:
        if directory[-1] != '/':
            directory += '/'
    """
    Reading the data
    """
    fluid_name = field
    field      =   np.fromfile(directory+fluid_name+quantity+'{:d}.dat'.format(it))
    field = field.reshape(p.nx*N,p.nz*N)
    field = field[::N,::N]
    return field
        

# -------------------------------------------------------------------  
# write AMR Grid
# -------------------------------------------------------------------  


def write_AMRgrid(mesh,params):

    print("writing AMR GRID")
    path_grid='amr_grid.inp'

    grid=open(path_grid,'w')

    grid.write('1 \n')              # iformat/ format number = 1
    grid.write('0 \n')              # Grid style (regular = 0)
    grid.write('101 \n')            # coordsystem: 100 < spherical < 200 
    grid.write('0 \n')              # gridinfo
    grid.write('1 \t 1 \t 1 \n')    # incl x, incl y, incl z
    # radius colatitude and azimuth
    grid.write(str(mesh.nx)+ '\t'+ str(mesh.ny)+'\t'+ str(mesh.nz)+'\n') 
    for i in range(mesh.nx+1):
        grid.write(str(mesh.xm[i])+'\t') 
    grid.write('\n')
    for i in range(mesh.ny+1):
        grid.write(str(mesh.ym[i])+'\t')
    grid.write('\n')
    # azimuth 
    # if azimuth grid goes between -pi and pi, add pi
    if ( np.abs(params.xmax-pi)/pi < 1e-3 ):
        print('forcing azimuth between 0 and 2pi')
        mesh.zm += pi
    for i in range(mesh.nz+1):
        grid.write(str(mesh.zm[i])+'\t') # added pi
    grid.write('\n')
    grid.close()

    # -------------------------------------------------------------------  
# writing out wavelength 
# -------------------------------------------------------------------  
def write_wavelength(M):
    Nw=M.Nw
    wmin=M.wmin
    wmax=M.wmax
    Pw = (wmax/wmin)**(1.0/(Nw-1))
    waves = np.zeros(Nw)
    waves[0] = wmin
    for i in list(range(1, Nw)):
        waves[i]=wmin*Pw**i
    path = 'wavelength_micron.inp'
    wave = open(path,'w')
    wave.write(str(Nw)+'\n')
    for i in list(range(Nw)):
        wave.write(str(waves[i])+'\n')
    wave.close()
    return waves


def write_stars(M,p=None,Plot=False):
    Nw=M.Nw
    wmin=M.wmin
    wmax=M.wmax
    Pw = (wmax/wmin)**(1.0/(Nw-1))
    
    waves = np.zeros(Nw)
    waves[0] = wmin
    for i in list(range(1, Nw)):
        waves[i]=wmin*Pw**i
    
    with open('stars.inp','w+') as f:
        f.write('2\n')
        f.write('1 %d\n\n'%(Nw))
        f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(M.Rstar,M.Mstar,0.0,0.0,0.0))
        for value in waves:
            f.write('%13.6e\n'%(value))
        f.write('\n%13.6e\n'%(-M.Tstar))
    
    
# -------------------------------------------------------------------  
# writing radmc3d.inp
# -------------------------------------------------------------------  

def write_radmc3dinp(incl_dust = 1,
                     incl_lines = 0,
                     nphot = 1000000,
                     nphot_scat = 1000000,
                     nphot_spec = 1000000,
                     nphot_mono = 1000000,
                     istar_sphere = 0,        ### stars treated as spheres
                     scattering_mode_max = 0,
                     tgas_eq_tdust = 1,
                     modified_random_walk = 0,
                     setthreads=2):

    print(( 'writing radmc3d.inp out'))

    RADMCINP = open('radmc3d.inp','w')
    inplines = ["incl_dust = "+str(int(incl_dust))+"\n",
                "incl_lines = "+str(int(incl_lines))+"\n",
                "nphot = "+str(int(nphot))+"\n",
                "nphot_scat = "+str(int(nphot_scat))+"\n",
                "nphot_spec = "+str(int(nphot_spec))+"\n",
                "nphot_mono = "+str(int(nphot_mono))+"\n",
                "istar_sphere = "+str(int(istar_sphere))+"\n",
                "scattering_mode_max = "+str(int(scattering_mode_max))+"\n",
                "tgas_eq_tdust = "+str(int(tgas_eq_tdust))+"\n",
                "modified_random_walk = "+str(int(modified_random_walk))+"\n",
                "setthreads="+str(int(setthreads))+"\n",
                "mc_scat_maxtauabs=5.d0\n",
                "rto_style = 1\n",   # 3 is binary output
                "rto_single = 1\n"]

    RADMCINP.writelines(inplines)
    RADMCINP.close()



# -------------------------------------------------------------------  
# writing dust and gas densities
# -------------------------------------------------------------------  
def rho(params,mesh, s,i, j, k, St,species):

    Sigma = 1.0*s
    R  = mesh.xm[i]                # radius in cm
    Th = mesh.ym[j]                # colatitude
    PHI = mesh.zm[k]
    R_cyl = np.sin(Th)*R 
    Z_cyl = np.cos(Th)*R

    h = params.aspectratio * (R_cyl/(R_Scaling*au))**(params.flaringindex)
    h *= np.sqrt(1.*params.alpha/(1.*params.alpha+1.*St))
    
    H = h*R_cyl
    rho = Sigma / (np.sqrt(2.*pi) * H*erf(np.cos(np.min(mesh.ym))/(np.sqrt(2)*h))) 
    exp_arg = -0.5 * ( np.cos(Th) / h )**2.0
    rho *= np.exp(exp_arg)
       
    return rho

# -------------------------------------------------------------------  
# writing dust densities
# -------------------------------------------------------------------  

def write_dust_densities(D, G, M, params, mesh,species, a, ndust):
    
    sigma_gas = G 
    sigma_dust = D
    ncells = mesh.nx*mesh.ny*mesh.nz
                
    # -------------------------------------------------------------------  
    # unit conversion
    # -------------------------------------------------------------------  
    if species==1:
        sigma_gas *= 2.0e33 / (au *au)
        sigma_gas *= 1.0 / (R_Scaling * R_Scaling)
    sigma_dust *= 2.0e33 / (au *au) 
    sigma_dust *= 1.0 / (R_Scaling * R_Scaling)
    if (CutEdges == True):
        Y_inf = params.ymin*params.dampingzone**(2.0/3.0)
        Y_sup = params.ymax*params.dampingzone**(2.0/3.0)
        Y_inf *= R_Scaling*au
        Y_sup *= R_Scaling*au
    # fractions
    print('writing dust and gas density fields')
    
    if (species == 1):
        DUSTOUT    = open('dust_density.inp','w')
    else:
        DUSTOUT    = open('dust_density.inp','a')
            
    if(species==1):
        DUSTOUT.write('1 \n')            # iformat  
        DUSTOUT.write(str(ncells)+' \n') # n cells
        DUSTOUT.write(str(ndust)+' \n')      # n species
            
        
    # DENSITIES 
    # only one dust species
    if ndust == 1:
        species = 2
    for k in range(mesh.nz):
        for j in range(mesh.ny):
            for i in range(mesh.nx):
                sigmag_ik = sigma_gas[i,k]
                sigmad_ik = sigma_dust[i,k]
                indo = (mesh.nx-N_outerdisk//N-20//N)
                #this is interpolating the density outside the hydro grid
                if i > indo:
                    sigmag_ik = sigma_gas[indo,0]*(mesh.xm[i]/mesh.xm[indo])**(-params.sigmaslope-0.5)
                    sigmad_ik = sigma_dust[indo,0]*(mesh.xm[i]/mesh.xm[indo])**(-params.sigmaslope-0.5)
                #Calculate local Stokes number needed for vertical analytical spreading
                if(sigmag_ik != 0.0):
                    St = pi/2.0*a*params.rho_solid/sigmag_ik
                else:
                    St= pi/2.0*a*params.rho_solid/(sigma_gas[200//N].mean(axis=-1)*(mesh.xm[i]/(0.4*R_Scaling*au))**(-0.5))
                #Call vertical spreading that translates 2d Sigma to 3d rho 
                rhod_ijk = rho(params,mesh,sigmad_ik,i,j,k,St,species)
                DUSTOUT.write(str(rhod_ijk)+' \n')
    if(species == ndust):
        DUSTOUT.close()
def write_dust_temp(temp,params,mesh,species,ndust):
    cs = temp
    ncells = mesh.nx*mesh.ny*mesh.nz
    if (species == 1):
        cs *= np.sqrt(G*2e33/(R_Scaling*au)) #convert to cgs 
        TOUT    = open('dust_temperature.dat','w')
    else:
        TOUT    = open('dust_temperature.dat','a')
    T = cs**2*2.3*m_H/kB
    print('Tmax=',T.max(),'Tmin=',T.min())
    if(species==1):
        TOUT.write('1 \n')            # iformat  
        TOUT.write(str(ncells)+' \n') # n cells
        TOUT.write(str(ndust)+' \n')      # n species

    for k in range(mesh.nz):
        for j in range(mesh.ny):
            for i in range(mesh.nx):
                TOUT.write(str(T[i,k])+' \n')
    if(species == ndust):
        TOUT.close()

def run_mctherm(modeldir=''):
    print(( "executing radmc3d mctherm \n"))
    os.system('radmc3d mctherm')
    
def run_raytracing(M,modeldir=''):
    command='radmc3d image lambda '+str(M.Lambda)+' npix '+str(M.npix)+' incl '+str(M.incl)+' posang '+str(M.posang)+' secondorder'

    #command='radmc3d image lambda '+str(M.Lambda)+' npix '+str(M.npix)+' incl '+str(M.incl)+' posang '+str(M.posang)+' sloppy '  # no secondorder is faster


    print(( command ))
    os.system(command)
        
def exportfits(M, Plot=False,Allwl=False):


#    outfile = 'image_i'+str(M.incl)+'_PA'+str(M.posang)+'.fits'
    outfile = 'image_out.fits'

    infile = 'image.out'

    # read header info:
    # iformat <=== For now this is 1 (or 2 for local observer mode)
    # im_nx im_ny
    # nlam
    # pixsize_x pixsize_y
    # lambda[1] ......... lambda[nlam+1]
    f = open(infile,'r')
    iformat = int(f.readline())
    im_nx, im_ny = tuple(np.array(f.readline().split(),dtype=int))
    nlam = int(f.readline()) 
    pixsize_x, pixsize_y = tuple(np.array(f.readline().split(),dtype=float))
    lbda = np.empty(nlam)
    for i in list(range(nlam)):
        lbda[i] = float(f.readline())
    f.readline()                # empty line

    # load image data
    images = np.loadtxt(infile, skiprows=(5+nlam))

    # calculate physical scales
    distance = M.distance    
    pixsize_x_deg = 180.0*pixsize_x / distance / pi
    pixsize_y_deg = 180.0*pixsize_y / distance / pi

    pixsurf_ster = pixsize_x*pixsize_x/distance/distance#pixsize_x_deg*pixsize_y_deg * (pi/180.)**2
    fluxfactor = 1e23 * pixsurf_ster
    
    if nlam>1:
        im = images.reshape(nlam,im_ny,im_nx)
    else:
        im = images.reshape(im_ny,im_nx)
        outfile = 'image_out_wl'+str(M.Lambda)+'.fits'

    if (Allwl):
        outfile = 'image_out_allwl.fits'
    
    if Plot:
        import matplotlib
        matplotlib.use('Agg')

        import matplotlib.pyplot as plt
        plt.imshow(im, cmap = 'plasma', origin='lower',aspect='auto')
        plt.axis('equal')
        plt.show()

    from astropy.io import fits

    hdu = fits.PrimaryHDU()
    hdu.header['BITPIX'] = '-32'
    hdu.header['NAXIS'] = 2
    hdu.header['NAXIS1'] = im_nx
    hdu.header['NAXIS2'] = im_ny
    hdu.header['BUNIT'] = 'JY/PIXEL'
    hdu.header['EPOCH'] = 2000.0
    hdu.header['EQUINOX'] = 2000.0
    hdu.header['LONPOLE'] = 180.0
    hdu.header['SPECSYS'] = 'LSRK    '
    hdu.header['CTYPE1'] = 'RA---TAN'
    hdu.header['CTYPE2'] = 'DEC--TAN'
    hdu.header['CRVAL1'] = 0.0
    hdu.header['CRVAL2'] = 0.0
    hdu.header['CDELT1'] = -1*pixsize_x_deg
    hdu.header['CDELT2'] = pixsize_y_deg
    hdu.header['CUNIT1'] = 'DEG     '
    hdu.header['CUNIT2'] = 'DEG     '
    hdu.header['CRPIX1'] = (im_nx+1.0)/2.0 
    hdu.header['CRPIX2'] = (im_ny+1.0)/2.0


    if nlam > 1:
  #      crpix3=np.floor((nlam-1)/2)+1
        crpix3=(nlam+1.0)/2.0
        restfreq = c * 1e4 / lbda[int(crpix3-1)] # micron to Hz
        nus = c * 1e4 / lbda                         # Hx 
        dvel = (lbda[1] - lbda[0])*c*1e-5/lbda[0]
        dnu = abs(nus[1] - nus[0])
        hdu.header['NAXIS'] = 3
        hdu.header['CTYPE3'] = 'FREQ    '
        hdu.header['CUNIT3'] = 'Hz      '
        hdu.header['CRPIX3'] = crpix3
        hdu.header['CRVAL3'] = restfreq
        hdu.header['CDELT3'] = dnu
        hdu.header['RESTFREQ'] = restfreq

    hdu.data = im*fluxfactor
    hdu.writeto(outfile, output_verify='fix', overwrite=True)
