import numpy as np
from astropy.io import fits
from astropy.convolution import convolve

def Gauss_filter(img, stdev_x, stdev_y, PA, Plot=False):
        ''' 
                img: image array
                stdev_x: float
                    BMAJ sigma dev.
                stdev_y: float
                    BMIN sigma dev.
                PA: East of North in degrees  '''

        image = img # fits.open(img)

        (nx0, ny0) = image.shape

        print( "nx ",nx0)
        print( "ny ",ny0)

        nx = min( nx0, int(7*stdev_x))
        nx = nx + ((nx+1) % 2)
        ny = nx
        print( "nx ",nx)

        x=np.arange(1,nx+1)
        print( "x shape ", x.shape)
        y=np.arange(1,ny+1)
        X, Y = np.meshgrid(x, y)
        print( "X shape ", X.shape)
        X0 = np.floor(nx/2)+1
        Y0 = np.floor(ny/2)+1
        X0 = nx/2
        Y0 = ny/2


        data = image

        #----------
        theta =   np.pi * (PA + 90.) / 180.
        A = 1
        a = np.cos(theta)**2/(2*stdev_x**2) + np.sin(theta)**2/(2*stdev_y**2)
        b = np.sin(2*theta)/(4*stdev_x**2) - np.sin(2*theta)/(4*stdev_y**2)
        c = np.sin(theta)**2/(2*stdev_x**2) + np.cos(theta)**2/(2*stdev_y**2)

        Z=A*np.exp(-(a*(X-X0)**2-2*b*(X-X0)*(Y-Y0)+c*(Y-Y0)**2))

        Z /= np.sum(Z)


        if Plot==True:
                import matplotlib.pyplot as plt
                import matplotlib.cm as cm

                plt.imshow(Z, cmap=cm.jet)
                plt.show()

        # Convolucion
        print( "data dims",data.shape)
        print( "kernel  dims",Z.shape)
        result = convolve(data, Z,boundary='extend')

        if Plot==True:
                plt.imshow(result, cmap=cm.magma)
                plt.show()

        return result
