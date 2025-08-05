import sys
import numpy as np
import os
import os.path
from scipy import ndimage
from astropy.io import fits as pf
import re
#include_path='/Users/simon/common/python/include/'
#sys.path.append(include_path)
from Resamp import *
from copy import deepcopy
from astropy.wcs import WCS

from pylab import *
#import matplotlib as plt

#import pyfits
#import matplotlib.pyplot as plt
#import re
#import scipy as S
#import scipy.ndimage
#from scipy import ndimage


def cartesian2polar(outcoords, inputshape, origin, fieldscale=1.):
    """Coordinate transform for converting a polar array to Cartesian coordinates. 
    inputshape is a tuple containing the shape of the polar array. origin is a
    tuple containing the x and y indices of where the origin should be in the
    output array."""

    rindex, thetaindex = outcoords
    x0, y0 = origin
    #print(origin)
    theta = thetaindex * 2 * np.pi / (inputshape[0]-1)
    #theta = 2. * np.pi - theta
    
    y = rindex*np.cos(theta)/fieldscale
    x = rindex*np.sin(theta)/fieldscale

    #print "r",rindex,"theta",theta,"x",x,"y",y
    
    ix = -x + x0
    iy = y +  y0

    return (iy,ix)


def polar2cartesian(outcoords, inputshape, origin, fieldscale=1.):
    yindex, xindex = outcoords
    x0, y0 = origin
    #print(origin)
#    theta0, r0 = origin
    nx=inputshape[0]-1
    ny=inputshape[1]-1
#    x0=((nx+1)/2)-1
#    y0=((ny+1)/2)-1
    x = -float(xindex) + x0
    y = float(yindex) -  y0
    
    #    theta = np.arctan2(-y, -x) + np.pi
    
    theta = np.arctan2(x, y) 
    if (theta < 0):
        theta = theta + 2.*np.pi
    #    theta = np.arctan2(-y, x) + np.pi
     
    thetaindex = (theta * nx / (2. * np.pi))
    #theta_index = np.round((theta + np.pi) * (inputshape[1]-1) / (2 * np.pi))
    
    #    thetaindex = np.arctan2(yindex-y0,xindex-x0) * nx / (2. * np.pi)
    rindex = fieldscale * np.sqrt(x**2 + y**2)
    
#    itheta = thetaindex + theta0
#    ir = rindex + r0

    #print xindex, yindex, " x ",x," y ", y," theta ",theta
    
    return (rindex,thetaindex)


def exec_polar_expansions(filename_source, workdir, PA, cosi, RA=False, DEC=False, alpha_min=False, Delta_min=False,
                          XCheckInv=False,DoRadialProfile=True,ProfileExtractRadius=-1,DoAzimuthalProfile=False,
                          PlotRadialProfile=True,a_min=-1,a_max=-1,zoomfactor=1.,Grid=False,y_label='',
                          ForceCube2Im=False,wBaseNoise=False,noise_radius=0.,wBaseNoiseCore=False):

    fieldscale=2. # shrink radial field of view of polar maps by this factor
    
    os.system("rm -rf mkdir "+workdir)
    
    os.system("mkdir "+workdir)
    inbasename=os.path.basename(filename_source)
    filename_fullim=re.sub('.fits', '_fullim.fits', inbasename)
    filename_fullim=workdir+filename_fullim
    
    hdu0 = pf.open(filename_source)
    hdr0 = hdu0[0].header
    if ((hdr0['NAXIS'] > 2) or ForceCube2Im):
        hdu=cube2im(filename_source,filename_fullim)
        im1=hdu.data
        hdr1=hdu.header
    else:
        os.system("rsync -va "+filename_source+" "+filename_fullim)
        hdu = pf.open(filename_fullim)
        im1=hdu[0].data
        hdr1=hdu[0].header

    hdr1.pop('CRVAL3', None)  
    

    if (not (isinstance(Delta_min,bool))):
        if (isinstance(RA,bool)):
            if (not RA):
                RA=hdr1['CRVAL1']
                DEC=hdr1['CRVAL2']
        
        RA=RA+(np.sin(alpha_min*np.pi/180.)*Delta_min/3600.)/np.cos(DEC*np.pi/180.)
        ##print "RA =", RA
        DEC=DEC+np.cos(alpha_min*np.pi/180.)*Delta_min/3600.
        ##print "DEC =",DEC

    elif (not isinstance(RA,float)):
        sys.exit("must provide a center")


    nx=int(hdr1['NAXIS1']/zoomfactor)
    ny=nx
    print(nx,ny)
    if ( (nx % 2) == 0):
        nx=nx+1
        ny=ny+1

    hdr2 = deepcopy(hdr1)

    

    hdr2['NAXIS1']=nx
    hdr2['NAXIS2']=ny
    hdr2['CRPIX1']=(nx+1)/2
    hdr2['CRPIX2']=(ny+1)/2
    hdr2['CRVAL1']=RA
    hdr2['CRVAL2']=DEC

    #print "filename_fullim ", filename_fullim
    
    resamp=gridding(filename_fullim,hdr2, fullWCS=False)


    fileout_centered=re.sub('fullim.fits', 'centered.fits', filename_fullim)
    pf.writeto(fileout_centered,resamp, hdr2, overwrite=True)

    #rotangel= PA - 180.
    rotangle= PA - 180.
    rotangle= PA
    im1rot = ndimage.rotate(resamp, rotangle, reshape=False)

    fileout_rotated=re.sub('fullim.fits', 'rotated.fits', filename_fullim)
    pf.writeto(fileout_rotated,im1rot, hdr2, overwrite=True)
    

    hdr3 = deepcopy(hdr2)
    hdr3['CDELT1']=hdr3['CDELT1']*cosi
    
    fileout_stretched=re.sub('fullim.fits', 'stretched.fits', filename_fullim)
    im3=gridding(fileout_rotated,hdr3)
    pf.writeto(fileout_stretched,im3, hdr2, overwrite=True)


    im_polar = sp.ndimage.geometric_transform(im3,cartesian2polar, 
                                              order=1,
                                              output_shape = (im3.shape[0], im3.shape[1]),
                                              extra_keywords = {'inputshape':im3.shape,'fieldscale':fieldscale,
                                                                'origin':(((float(nx-1))/2.0),((float(ny-1))/2.0))})
                                                                #'origin':(((nx+1)/2)-1,((ny+1)/2)-1)}) 
    
    nphis,nrs=im_polar.shape
    
    hdupolar = pf.PrimaryHDU()
    hdupolar.data = im_polar
    hdrpolar=hdupolar.header
    hdrpolar['CRPIX1']=1
    hdrpolar['CRVAL1']=0.
    hdrpolar['CDELT1']=2. * np.pi / nphis
    hdrpolar['CRPIX2']=1
    hdrpolar['CRVAL2']=0.
    hdrpolar['CDELT2']=(hdr3['CDELT2'] / fieldscale)
    hdupolar.header = hdrpolar
    
    fileout_polar=re.sub('fullim.fits', 'polar.fits', filename_fullim)
    hdupolar.writeto(fileout_polar, overwrite=True)

    ######################################################################
    # profiles

    if (DoRadialProfile):

        Iprof = np.average(im_polar,axis=1)
        sIprof = np.std(im_polar,axis=1)
        rrs =  3600.*(np.arange(hdrpolar['NAXIS2'])-hdrpolar['CRPIX2']+1)*hdrpolar['CDELT2']+hdrpolar['CRVAL2']
        #np.savetxt('test.out', (rrs,Iprof))  

        bmaj=hdr2['CDELT2']
        if ('BMAJ' in hdr2):
            bmaj=hdr2['BMAJ']
        ##print "bmaj = ",bmaj,"\n";
        Nind=2.*np.pi*rrs * cosi  /(bmaj*3600.) #cosi *
        iNind1=np.argmin(np.fabs(Nind - 1.0))
        rNind1=rrs[iNind1]
        ##print "radius for Nind=1 ", rNind1
        ##print "rrs ",
        Nind[np.where(Nind < 1.0)]=1.0
        dispIprof=sIprof.copy()

        if (noise_radius <= 0.):
            noise_radius=rNind1
        
        if (wBaseNoise):
            noise_basal=dispIprof[np.argmin(np.fabs(rrs - noise_radius))]
            if (wBaseNoiseCore):
                #dispIprof[np.where( (rrs < noise_radius) & (dispIprof<noise_basal))] = noise_basal
                dispIprof[np.where(rrs < noise_radius)] = noise_basal
            else:
                dispIprof[np.where(dispIprof<noise_basal)] = noise_basal
            ##print ">using noise floor of ",noise_basal,"  from dispersion at radius ", noise_radius
        else:
            print( ">not applying noise floor")
            
        sIprof = dispIprof/np.sqrt(Nind)
        
        save_prof = np.zeros((hdrpolar['NAXIS2'],4))
        #print save_prof.shape
        save_prof[:,0] = rrs
        save_prof[:,1] = Iprof
        save_prof[:,2] = sIprof
        save_prof[:,3] = dispIprof
        fileout_radialprofile=re.sub('fullim.fits', 'radial_profile.dat', filename_fullim)
        np.savetxt(fileout_radialprofile, save_prof)   # x,y,z equal sized 1D arrays


        ######################################################################
        # Polar average
        
    
        Iprof_b = Iprof.reshape(Iprof.shape[0],1)    
        #im_polar_av = im_polar.copy
        im_polar_av = np.copy(im_polar)
        
        im_polar_av[:,:] = Iprof_b
        
        
        fileout_polar_av=re.sub('fullim.fits', 'polar_av.fits', filename_fullim)
        hdupolar.data = im_polar_av
        hdupolar.writeto(fileout_polar_av, overwrite=True)
    
    
        ######################################################################
        # AZIM IM 
        
        
        imazim = sp.ndimage.geometric_transform(im_polar_av,polar2cartesian,
                                                order=0,
                                                output_shape = (im_polar.shape[0], im_polar.shape[1]),
                                                extra_keywords = {'inputshape':im_polar.shape,'fieldscale':fieldscale,
                                                                  'origin':(((float(nx-1))/2.0),((float(ny-1))/2.0))})
                                                                  #'origin':(((nx+1)/2)-1,((ny+1)/2)-1)}) 
        
    
        fileout_stretched_av=re.sub('fullim.fits', 'stretched_av.fits', filename_fullim)
        pf.writeto(fileout_stretched_av,imazim, hdr2, overwrite=True)


        ######################################################################
        # back to sky - project 

    
        hdr3 = deepcopy(hdr2)
        hdr3['CDELT1']=hdr3['CDELT1']/cosi
        
        fileout_proj=re.sub('fullim.fits', 'azim_av_proj.fits', filename_fullim)
        im4=gridding(fileout_stretched_av,hdr3)
        pf.writeto(fileout_proj,im4, hdr2, overwrite=True)
        
    
        ######################################################################
        # back to sky - rotate
        
        im4drot = ndimage.rotate(im4, -rotangle, reshape=False)
        
        fileout_drotated=re.sub('fullim.fits', 'azim_av_drot.fits', filename_fullim)
        pf.writeto(fileout_drotated,im4drot, hdr2, overwrite=True)
        
    
        ######################################################################
        # azimuthal profiles

    if ( (DoAzimuthalProfile) or (ProfileExtractRadius  > 0)):


        
        im_polar_4profiles = np.copy(im_polar)

        rs =  3600.*((np.arange(hdrpolar['NAXIS2'])-hdrpolar['CRPIX2']+1)*hdrpolar['CDELT2']+hdrpolar['CRVAL2'])
        if (a_min > 0):
            ia_min = np.argmin(np.abs(rs - a_min))
            ##print "a_min",a_min," ia_min",ia_min
            im_polar_4profiles[0:ia_min,:] = 0

        if (a_max > 0):
            ia_max = np.argmin(np.abs(rs - a_max))
            ##print "a_max",a_max," ia_min",ia_max
            im_polar_4profiles[ia_max:,:] = 0

        #im_polar_4profiles[0:20,:] = 0

        #import matplotlib.pyplot as plt
        #plt.imshow(im_polar_4profiles)
        #plt.show()

        Imax_profile=np.amax(im_polar_4profiles,axis=0)
        phis =  (180./np.pi)*((np.arange(hdrpolar['NAXIS1'])-hdrpolar['CRPIX1']+1)*hdrpolar['CDELT1']+hdrpolar['CRVAL1'])

        (nphis, nrs) = im_polar_4profiles.shape

        #print "nphis ",nphis,"len phis",len(phis)
        #print "nrs ",nrs
        #iphis=np.arange(1,nphis+1)
        #irs=np.arange(1,nrs+1)
        #iiphis,iirs  = np.meshgrid(iphis, irs)
        #rrs =  3600.*(iirs-hdrpolar['CRPIX2']+1)*hdrpolar['CDELT2']+hdrpolar['CRVAL2']
        #pphis =  (180./np.pi)*(iiphis-hdrpolar['CRPIX1']+1)*hdrpolar['CDELT1']+hdrpolar['CRVAL1']


        ivec_rmax=np.argmax(im_polar_4profiles,axis=0)
        Rmax_profile =  3600.*( (ivec_rmax-hdrpolar['CRPIX2']+1)*hdrpolar['CDELT2']+hdrpolar['CRVAL2'])
        ##print "Rmax shape:",Rmax_profile.shape
        
        
        
        WholeRing=False
        if (ProfileExtractRadius > 0):
            ringrad=ProfileExtractRadius
        elif (WholeRing):
            #Iprof = np.average(im_polar_4profiles,axis=1)
            #ringrad=np.asscalar(Iprof[np.argmax(Iprof)])
            ringrad=np.sum(Rmax_profile*Imax_profile)/np.sum(Imax_profile)
            ##print "Rcav ", ringrad
        else:
            iphimax = np.argmax(Imax_profile)
            ##print "iphimax",iphimax
            iringrad=ivec_rmax[iphimax]
            ringrad=Rmax_profile[iphimax]

        ##print "iringrad ",iringrad," ringrad ",ringrad 

        Icyl_profile=im_polar_4profiles[iringrad,:]

        
        
        #np.savetxt(workdir+'rmax_profile.dat', (rs,Rmax_profile))  
        #np.savetxt(workdir+'profile_peak.dat', (rs,Imax_profile))  
        #np.savetxt(workdir+'profile_cyl.dat', (rs,Icyl_profile))  
        
        save_prof = np.zeros((hdrpolar['NAXIS2'],2))
        save_prof[:,0] = phis
        save_prof[:,1] = Rmax_profile
        np.savetxt(workdir+'rmax_profile.dat', save_prof)  
        
        save_prof = np.zeros((hdrpolar['NAXIS2'],2))
        save_prof[:,0] = phis
        save_prof[:,1] = Imax_profile
        np.savetxt(workdir+'profile_peak.dat', save_prof)  

        save_prof = np.zeros((hdrpolar['NAXIS2'],2))
        save_prof[:,0] = phis
        save_prof[:,1] = Icyl_profile
        np.savetxt(workdir+'profile_cyl.dat', save_prof)  

        avring=np.sum(Rmax_profile*Imax_profile)/np.sum(Imax_profile)
        #rmsring=np.std(Rmax_profile*Imax_profile/np.sum(Imax_profile))
        rmsring2=np.sqrt(np.sum( (Rmax_profile-avring)**2*Imax_profile/np.sum(Imax_profile)))
        #print "rmsring ",rmsring,"rmsring2 ", rmsring2 
        if (Grid):
            return rmsring2 
        
        

    
    
    ######################################################################
    # CROSS CHECK INVERSE TRANSFORM

    if (XCheckInv):
        im_x = sp.ndimage.geometric_transform(im_polar,polar2cartesian,
                                              order=0,
                                              output_shape = (im_polar.shape[0], im_polar.shape[1]),
                                              extra_keywords = {'inputshape':im_polar.shape,'fieldscale':fieldscale,
                                                                'origin':((float(nx-1)/2),(float(ny-1)/2))}) 
        
        fileout_stretched_x=re.sub('fullim.fits', 'stretched_x.fits', filename_fullim)
        pf.writeto(fileout_stretched_x,im_x, hdr2, overwrite=True)
        
        hdr3 = deepcopy(hdr2)
        hdr3['CDELT1']=hdr3['CDELT1']/cosi
        fileout_proj_x=re.sub('fullim.fits', 'x_proj.fits', filename_fullim)
        im4_x=gridding(fileout_stretched_x,hdr3)
        pf.writeto(fileout_proj_x,im4_x, hdr2, overwrite=True)
        
        im4_x_drot = ndimage.rotate(im4_x, -rotangle, reshape=False)
        fileout_drotated_x=re.sub('fullim.fits', 'x_drot.fits', filename_fullim)
        pf.writeto(fileout_drotated_x,im4_x_drot, hdr2, overwrite=True)

        
        fileout_diff_x=re.sub('fullim.fits', 'x_diff.fits', filename_fullim)
        pf.writeto(fileout_diff_x,resamp-im4_x_drot, hdr2, overwrite=True)

    
    
    
    ######################################################################
    # profile plotting
    import matplotlib.pyplot as plt

    if (PlotRadialProfile):
        #import matplotlib as plt
        #### from pylab import *
        import matplotlib.colors as colors
        
        # -----------------------------------------------------------
        # nice fonts
        # -----------------------------------------------------------
        matplotlib.rc('font', family='sans-serif') 
        matplotlib.rcParams.update({'font.size': 12})

        # normalization
        sIprof = sIprof/np.max(Iprof)
        dispIprof = dispIprof/np.max(Iprof)
        Iprof = Iprof/np.max(Iprof)
        
        fig = plt.figure()
        #figsize(7, 5)
        axprofile = fig.add_subplot(111)
        new_rrs = rrs
        rmax=np.max(rrs)
        axprofile.get_xticklabels()
        axprofile.get_yticklabels()
        plt.xlim(0.,0.9)
        plt.ylim(-0.1,1.1)
        plt.axvspan(0, 0.1, facecolor='gray', alpha=0.4)
        plt.plot(new_rrs,new_rrs*0.,color='black',linewidth=0.1,linestyle='solid')
        #plt.plot(rrs,Icutav,color='blue',linewidth=1,linestyle='solid')
        plt.plot(new_rrs,Iprof,color='grey',linewidth=0.1,linestyle='solid')
        plt.fill_between(new_rrs, Iprof+sIprof, Iprof-sIprof, lw=0.1,color='r', alpha=0.3, interpolate=True, step='mid')
        plt.fill_between(new_rrs, Iprof+dispIprof, Iprof-dispIprof, lw=0.1,color='b', alpha=0.3, interpolate=True)
        plt.ylabel('Normalized Surface Brightness')
        plt.xlabel(r'$r$ / arcsec')
        
        fileout_fig=re.sub('fullim.fits', 'fig_profile.pdf', filename_fullim)
        plt.savefig(fileout_fig, bbox_inches='tight')
   
        new_rrs_au = (np.pi / 180.) * (rrs /3600.) * 14933571.969912
        return [Iprof, sIprof, dispIprof, new_rrs, new_rrs_au]
