#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 14:15:20 2018

@author: Henrik Åhl (henrik.aahl@slcu.cam.ac.uk)

-------------------------------------------------------------------------------

This code is produced for the publication "Quantitative analysis of auxin sensing 
in leaf primordia argues against proposed role in regulating leaf dorsoventrality",
submitted to eLife 18-06-2018, and authored by:

- Neha Bhatia (bhatia@mpipz.mpg.de)
- Henrik Åhl (henrik.aahl@slcu.cam.ac.uk)
- Henrik Jönsson (henrik.jonsson@slcu.cam.ac.uk)
- Marcus Heisler (marcus.heisler@sydney.edu.au)

For queries relating to the paper, contact Marcus Heisler (marcus.heisler@sydney.edu.au).
Questions related to the code are best addressed to Henrik Åhl (henrik.aahl@slcu.cam.ac.uk).

-------------------------------------------------------------------------------

This particular code file deals with the deconvolution of the image files. Note
that the files are large numpy arrays, which consume a lot of memory when processed, 
and can cause issues on smaller desktop computers.

By default, the script takes data from <path-to-repo>/intermediate_data/cropped_leaves, 
and saves the output in <path-to-repo>/intermediate_data/cropped_and_deconvolved_leaves.
The laser settings used for the deconvolution are stored in <path-to-repo>/misc/laser_settings.txt.

"""


import os
project_path = '/home/henrik/projects/r2d2_leaves/'
os.chdir(os.path.join(project_path, 'code'))

import numpy as np
import pandas as pd
import tifffile as tiff
from skimage.restoration import richardson_lucy
from scipy.signal import wiener
from multiprocessing import Pool
import gc

from misc import listdir, autocrop, get_resolution
from external import psf

# Get paths
data_dir = os.path.join(project_path, 'intermediate_data', 'cropped_leaves')
files = listdir(data_dir, include='.tif')
files.sort()

# Set parameters
MAGNIFICATION = 25
NCHANNELS = 3
NA = 0.95
DECONVOLUTION_ITERATIONS = [7, 8, 8]
NPROC = 20

# Get laser settings used
laser_settings = pd.read_csv(os.path.join(project_path, 'misc', 'laser_settings.txt'), 
                             sep='\t')

def create_psf(zshape, rshape, zdims, rdims, ex_wavelen, em_wavelen, pinhole_radius,
               num_aperture=1.0, refr_index=1.333,  pinhole_shape='square',
               psf_type=(psf.ISOTROPIC | psf.CONFOCAL), magnification=20.0):
    """ 
    Create a point-spread function (PSF) for using the Richards and Wolf formula.
    
    See psf.PSF function for a detailed description.
    
    Parameters
    ----------
    zshape : int
        Extent of the Z-dimension of the image of interest.
    rshape : int
        Extent of the radial dimension of the image of interest.
    zdims : float
        Physical Z-extent. Value given in um.
    rdims : float
        Physical R-extent. Value given in um.
    ex_wavelen : float
        Excitation wavelength of corresponding fluorophore. Value given in nm.
    em_wavelen : float
        Emission wavelength of corresponding fluorophore. Typically taken as the
        mean of the emission range. Value given in nm.
    pinhole_radius : float
        Radius of the microscope pinhole. Value given in um.
    num_aperture : float
        Numerical aperture of the microscope lens.
    refr_index : float
        Refractive index of the imaging solution.
    pinhole_shape : str
        One of 'square' or 'round',
    psf_type : int
        See psf.PSF. Default is (psf.ISOTROPIC | psf.CONFOCAL)
    magnification : float
        Magnification factor of the microscope lens.
    
    Returns
    -------
    obsvol : 3D np.ndarray
        Image array depicting the calculated PSF. 
    
    """

    args = dict(shape=(int(zshape), int(rshape)),
                dims=(zdims, rdims),
                ex_wavelen=ex_wavelen,
                em_wavelen=em_wavelen,
                num_aperture=num_aperture,
                refr_index=refr_index,
                pinhole_radius=pinhole_radius,
                pinhole_shape=pinhole_shape,
                magnification=magnification)

    obsvol = psf.PSF(psf_type, **args).volume()

    return obsvol

def deconvolve(data, psf_vol, iterations, threshold=0, wiener_size=(3, 3, 3), wiener_noise=500000):
    """ 
    Deconvolve an image using the Richardson Lucy iterative algorithm. Subsequently
    apply a Wiener filter for denoising and improving the signal further.
    
    data : np.array
        Input image.
    iterations : int
        Number of iterations to deconvolve for.
    threshold : float
        Background filtering threshold. Intensities with a lower value are set to 0.
    wiener_size : np.array 
        Wiener filter footprint of same dimensionality as `data`.
    wiener_noise : float
        Wiener filter noise parameter.
        
    Returns
    -------
    data : np.array
        The deconvolved image.
    
    """
    
    data[data < threshold] = 0
    data = richardson_lucy(image=data, psf=psf_vol, iterations=iterations, clip=False)
    gc.collect()

    data[data < 0] = 0
    if wiener_noise != 0:
        data = wiener(data, mysize=wiener_size, noise=wiener_noise)
    data[data < 0] = 0

    return data

#for fname in files:
def deconvolve_single(fname):
    """ Convenience function for deconvolving a single file. """
    
    # Get data and relevant parameters
    sample, leaf, side = os.path.basename(fname).split('_')[0:3]
    sample = int(sample)
    f = tiff.TiffFile(fname)
    data = f.asarray()
    data = autocrop(data, fct=np.max, threshold=0)
    data = data.astype(np.float64)
    resolution = get_resolution(f)
    del f
    
    # Get dimensions
    zshape = data.shape[0] // 2. + 2.
    rshape =  data.shape[-1] // 2. + 2.
    zdims = zshape * resolution[0]
    rdims = rshape * resolution[1]

    # Get laser values
    sampledata = laser_settings.loc[laser_settings['sample'] == sample]
    ex_wavelens = np.array([sampledata['excitation_wavelength_1'].values[0],     
                            sampledata['excitation_wavelength_2'].values[0], 
                            sampledata['excitation_wavelength_3'].values[0]])
    em_wavelens = np.array([(float(sampledata['range_low_nm_1'].values[0]) + 
                                 sampledata['range_high_nm_1'].values[0]) / 2,
                            (float(sampledata['range_low_nm_2'].values[0]) + 
                                 sampledata['range_high_nm_2'].values[0]) / 2,
                            (float(sampledata['range_low_nm_3'].values[0]) + 
                                 sampledata['range_high_nm_3'].values[0]) / 2])
    pinhole_radii = np.array([sampledata['pinhole_um'].values[0] / 2 / MAGNIFICATION,
                              sampledata['pinhole_um'].values[0] / 2 / MAGNIFICATION,
                              sampledata['pinhole_um'].values[0] / 2 / MAGNIFICATION])

    # Deconvolve
    for ii in xrange(NCHANNELS):
        gc.collect()
        cpsf = create_psf(zshape, rshape, zdims, rdims, ex_wavelen=ex_wavelens[ii],
                          em_wavelen=em_wavelens[ii], pinhole_radius=pinhole_radii[ii],
                          num_aperture=NA, magnification=MAGNIFICATION, pinhole_shape='round')
        data[:, ii] = deconvolve(data[:, ii], psf_vol=cpsf, iterations=DECONVOLUTION_ITERATIONS[ii], threshold=150)
    data[data > np.iinfo(np.uint16).max] = np.iinfo(np.uint16).max
    data[data < 0] = 0

    # Save
    tiff.imsave(os.path.join(project_path, 'intermediate_data', 'cropped_and_deconvolved_leaves', 
                             os.path.splitext(os.path.basename(fname))[0] + '_deconvolved.tif'), 
        data=data.astype(np.uint16), shape=data.shape, dtype=data.dtype, imagej=True, 
        metadata={'spacing' : resolution[0]}, resolution=1. / resolution[1:])

# Run stuff
p = Pool(NPROC)
p.map(deconvolve_single, files)
p.close()
p.join()
