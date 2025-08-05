"""
Code interpolates a FITS image into the grid defined by another FITS header, using ndimage.map_coordinates from scipy.
"""
import sys
import scipy as sp
from scipy.ndimage import map_coordinates
from astropy.io import fits as pf
from astropy.wcs import WCS


def cube2im(filename,fileout=False):
    dcube1, hdr1 = datafits(filename)
    hdr1.pop('PC3_1', None)  
    hdr1.pop('PC4_1', None)  
    hdr1.pop('PC3_2', None)  
    hdr1.pop('PC4_2', None)  
    hdr1.pop('PC1_3', None)  
    hdr1.pop('PC2_3', None)  
    hdr1.pop('PC3_3', None)  
    hdr1.pop('PC4_3', None)  
    hdr1.pop('PC1_4', None)  
    hdr1.pop('PC2_4', None)  
    hdr1.pop('PC3_4', None)  
    hdr1.pop('PC4_4', None)
    hdr1.pop('PC03_01', None)
    hdr1.pop('PC04_01', None)
    hdr1.pop('PC03_02', None)
    hdr1.pop('PC04_02', None)
    hdr1.pop('PC01_03', None)
    hdr1.pop('PC02_03', None)
    hdr1.pop('PC03_03', None)
    hdr1.pop('PC04_03', None)
    hdr1.pop('PC01_04', None)
    hdr1.pop('PC02_04', None)
    hdr1.pop('PC03_04', None)
    hdr1.pop('PC04_04', None)
    hdr1.pop('CTYPE3', None) 
    hdr1.pop('CRVAL3', None) 
    hdr1.pop('CDELT3', None) 
    hdr1.pop('CRPIX3', None) 
    hdr1.pop('CUNIT3', None) 
    hdr1.pop('CTYPE4', None) 
    hdr1.pop('CRVAL4', None) 
    hdr1.pop('CDELT4', None) 
    hdr1.pop('CRPIX4', None) 
    hdr1.pop('CUNIT4', None) 
    hdr1.pop('OBJECT', None)
    hdr1.pop('PC01_01', None)
    hdr1.pop('PC02_01', None)
    hdr1.pop('PC01_02', None)
    hdr1.pop('PC02_02', None)
    hdr1.pop('HISTORY', None)
    hdr1.pop('COMMENT', None)
    
    hdr1.pop('CROTA3', None)
    hdr1.pop('CROTA4', None)

    hdr1.pop('', None) 
    hdr1.pop('TELESCOP', None) 
    hdr1.pop('LONPOLE', None) 
    hdr1.pop('LATPOLE', None) 
 

    if (len(dcube1.shape) > 3):
        im1=dcube1[0,0,:,:]
    elif (len(dcube1.shape) > 2):
        im1=dcube1[0,:,:]
    else:
        im1=dcube1


        
    im1=sp.nan_to_num(im1)

    if (isinstance(fileout,str)):
        pf.writeto(fileout,im1, hdr1, clobber=True)

    hdu = pf.PrimaryHDU()
    hdu.data = im1
    hdu.header = hdr1
    return hdu


def datafits(namefile):
    """
    Open a FITS image and return datacube and header.
    """
    datacube = pf.open(namefile)[0].data
    hdr = pf.open(namefile)[0].header
    return datacube, hdr


def gridding(imagefile_1, imagefile_2,fileout=False,fullWCS=True):
    """
    Interpolates Using ndimage and astropy.wcs for coordinate system.
    """

    if (isinstance(imagefile_1,str)):
        im1, hdr1 = datafits(imagefile_1)
    elif (isinstance(imagefile_1,pf.hdu.image.PrimaryHDU)):
        im1 = imagefile_1.data
        hdr1 = imagefile_1.header
    elif (isinstance(imagefile_1,pf.hdu.hdulist.HDUList)):
        im1 = imagefile_1[0].data
        hdr1 = imagefile_1[0].header
    else:
        sys.exit("not an recognized input format")
        
    
    if (isinstance(imagefile_2,str)):
        im2, hdr2 = datafits(imagefile_2)
    else:
        hdr2=imagefile_2

    hdr1.pop('CRVAL3', None)  
    hdr2.pop('CRVAL3', None)  

        
    w1 = WCS(hdr1)
    w2 = WCS(hdr2)
    
    n2x = hdr2['NAXIS1']
    n2y = hdr2['NAXIS2']
    k2s=sp.arange(0,n2x)
    l2s=sp.arange(0,n2y)
    kk2s, ll2s = sp.meshgrid(k2s, l2s)

    if (fullWCS):
        xxs2wcs, yys2wcs = w2.all_pix2world(kk2s, ll2s, 0)
        kk1s, ll1s = w1.all_world2pix(xxs2wcs,yys2wcs,0,tolerance=1e-12)
    else:
        xxs2wcs, yys2wcs = w2.wcs_pix2world(kk2s, ll2s, 0)
        kk1s, ll1s = w1.wcs_world2pix(xxs2wcs,yys2wcs,0)


    im1=sp.nan_to_num(im1)
  
    resamp = map_coordinates(im1, [ll1s, kk1s],prefilter=False) #,order=1

    resamp=sp.nan_to_num(resamp)

    if (fileout):
        pf.writeto(fileout,resamp, hdr2, clobber=True)

    return resamp


    


