########################################################################
# A colletion of useful functions for reduction of VLT Sinfoini data.  #
#                                                                      #
#                    author: Ryleigh Davis                             #
########################################################################

from fits_dataclass import Image, ImageSet, Transform
import numpy as np
from photutils import centroid_2dg
import copy

### BASIC IMAGE UTILS ###

def get_wave(img):
    """ Use the image header values (from the official VLT Sinfoini data
    reduction pipeline) to identify the wavelength of each layer in the
    spectral cube.
    
    INPUT: 
        img -> Image: Image(data, header)
    
    OUTPUT: 
        wave -> arr: array containing the wavelength along the spectral 
                cube data contained in img
    """
    #TO DO: Add wavelength to Image dataclass i.e. spectral cube
    cen_lambda = img.header['CRVAL3']
    cen_pix = int(img.header['CRPIX3'])
    micron_per_pix = img.header['CDELT3']
    
    wave = np.zeros(img.shape[0])
    for i in range(img.shape[0]):
        wave[i] = cen_lambda + (i-cen_pix)*micron_per_pix
    
    return wave




### TELLURIC CORRECTION (VIA REFERENCE STAR) ###

def get_center(img):
    """ Return the center x,y coordinate of a reference star in 
    an image using gaussian centering."""

    x, y = centroid_2dg(img)
    return (x,y)

def telluric_correct(img, std, box_size=30):
    '''Perform telluric star division on an image.
    
    INPUTS:
        img: Image: image to be corrected (Europa)
        std: Image: standard star image
        box_size: int: size of box to sum std star'''
    
    #TO DO: Be able to scale absorption lines
    #        - separate continuum and absorption line removal
    
    #Get center of std star at mid-wavelength
    cen_wave = int(img.data.shape[0]/2)
    x,y=get_center(std.data[cen_wave])
    
    #Sum std star in a box of width box_size around the center
    summed = np.array([np.sum(std.data[i][int(x-box_size/2):int(x+box_size/2),
                                        int(y-box_size/2):int(y+box_size/2)]
                             ) for i in range(std.shape[0])])
    
    #Divide each spaxel by summed
    
    #TO DO: assert len(summed) == len(spaxel)
    
    
    div = np.zeros(img.data.shape)

    for i in range(img.data.shape[1]):
        for j in range(img.data.shape[2]):
        
            div[:,i,j] = img.data[:,i,j].data/summed
    
    #Add std star file name to img header
    h = copy.copy(img.header)
    h['StdCorr'] = std.header['ARCFILE']
         
    return Image(_data=div, header=h)



###  Combine corner and center images into a single image cube ###

def get_mask(imgs):
    
    med = np.nanmedian(np.dstack(imgs[0].data), axis=2)
    
    rolled = med - np.roll(med,2)
    
    mask = abs(rolled) < 0.4*np.std(rolled)

    #Mask out Border
    #TO DO: Maybe automate finding this region?
    mask[0:7,:] = False
    mask[:,0:4] = False
    mask[59:,:] = False
    mask[:,59:] = False
    
    return mask

def get_img_shift(RA1, DEC1, RA2, DEC2, PA=0, 
                  plate_scale=(0.0125, 0.0125)):
    """ Return x,y shift in pixel space for image based on RA
    and DEC pointing information."""
    
    if PA == 0:
        #Make sure PA is 0
    
        delta_RA = RA1 - RA2
        delta_DEC = DEC2 - DEC1
    
        x_shift = delta_RA*3600 #in arcsec
        y_shift = delta_DEC*3600 #in arcsec
        
        return x_shift/plate_scale[0], y_shift/plate_scale[1]
    
    else:
        #TO DO: update to handle different PA
        return None
    
def get_img_loc(x_shift, y_shift):
    """Return lower x and y coordinate for image on combined image"""
    
    return x_shift+32, y_shift+32

def get_comb_mask(imgs):
    """Get a combined mask
        0: no data
        1-4: corner data, img i
        5: center data
        6-9: corner and center data, 5+i"""

    comb_mask = np.zeros((128,128)) #Blank combined mask

    mask = get_mask(imgs)

    for i in range(1,len(imgs)):
        x_shift, y_shift = get_img_shift(imgs[0].header['CRVAL1'], 
                                      imgs[0].header['CRVAL2'], 
                                      imgs[i].header['CRVAL1'], 
                                      imgs[i].header['CRVAL2'])
    
        #TO DO: Update to be able to have subpixel accuracy
        x_shift, y_shift = int(x_shift), int(y_shift) 
        
        lx,ly = get_img_loc(x_shift, y_shift)
    
        comb_mask[ly:ly+64,lx:lx+64] = i*mask
     
    #Add in center img
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):

            if mask[i,j] == 1:
                if comb_mask[32+i,32+j] == 0:
                    comb_mask[32+i,32+j] = 5 #Take center data
                else:
                    comb_mask[32+i,32+j] = 5+comb_mask[32+i,32+j] #Corner and center data
            #otherwise take corner data (1-4 depending on img)
    
    return comb_mask
