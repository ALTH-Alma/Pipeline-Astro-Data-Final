#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import numpy as np
import math 
#import copy
import re
import sys
from Gausssmooth import *


#filename_ALMA='DoAr44_B7_zoom.fits'

filename_RTALMA='image_out_wl1300.fits'
filename_smooth='smooth_image_out_wl1300_smooth.fits'

#f = fits.open(filename_ALMA)
#cube = f[0].data
#hdr= f[0].header

fRT = fits.open(filename_RTALMA)
im = fRT[0].data
hdrRT= fRT[0].header

spike=np.zeros(im.shape)

spike[int(hdrRT['CRPIX2']-1), int(hdrRT['CRPIX1']-1)] = 1


###################################################################### 

bmaj = 1.633395751317E-05
bmin = 1.142132199473E-05      # clean_wrong
bpa = -8.714191436768E+01

bmaj = 0.0775101/3600 /2
bmin = 0.0514892/3600 /2    # clean_real
bpa = -85.936

bmaj = 1.714785686798E-05  /np.sqrt(3)
bmin = 1.532424034344E-05  /np.sqrt(3)    # uvmem_nogrid
bpa =  8.752417755127E+01  

######################################################################
#

# smooth to  clean beam

stdev_x = (bmaj/(2*np.sqrt(2*np.log(2)))) / hdrRT['CDELT2']
stdev_y = (bmin/(2*np.sqrt(2*np.log(2)))) / hdrRT['CDELT2']
PA=bpa

#print "stdev_x",stdev_x
#print "stdev_y",stdev_y
#print "PA",PA

im_smooth = Gauss_filter(im,stdev_x,stdev_y,PA,Plot=False)


beam = (np.pi/(4.*np.log(2))) * bmaj * bmin
pixel = hdrRT['CDELT2']**2


im_smooth *=   beam / pixel

hdu = fits.PrimaryHDU()
hdu.data = im_smooth
hdu.header = hdrRT

hdu.writeto(filename_smooth,overwrite=True)

