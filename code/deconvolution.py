#!/usr/bin/env python2
# -*- coding: utf-8 -*-
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
data_dir = os.path.join(project_path, 'intermediate_data/cropped_leaves')
files = listdir(data_dir, include='.tif')
files.sort()

# Set parameters
MAGNIFICATION = 25
NCHANNELS = 3
NA = 0.95
DECONVOLUTION_ITERATIONS = [7, 8, 8]
NPROC = 20

# Get laser settings used
laser_settings = pd.read_csv(os.path.join(project_path, 'misc/laser_settings.txt'), sep='\t')

def create_psf(zshape, rshape, zdims, rdims, ex_wavelen, em_wavelen, pinhole_radius,
               num_aperture=1.0, refr_index=1.333,  pinhole_shape='square',
               psf_type=(psf.ISOTROPIC | psf.CONFOCAL), magnification=20.0):

    args = dict(shape=(int(zshape), int(rshape)),
                dims=(zdims, rdims),
                ex_wavelen=ex_wavelen,
                em_wavelen=em_wavelen,
                num_aperture=num_aperture,
                refr_index=refr_index,
                pinhole_radius=pinhole_radius,
                pinhole_shape=pinhole_shape,
                magnification=magnification)

    obsvol = psf.PSF(psf_type, **args)

    return obsvol.volume()

def deconvolve(data, psf_vol, iterations, threshold=0):
    data[data < threshold] = 0
    data = richardson_lucy(image=data, psf=psf_vol, iterations=iterations, clip=False)
    gc.collect()

    data[data < 0] = 0
    data = wiener(data, mysize=(3, 3, 3), noise=500000)
    data[data < 0] = 0

    return data

#for fname in files:
def deconvolve_single(fname):
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
