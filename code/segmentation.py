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

This particular script deals with peforming the nuclear segmentation used for the 
quantitative analysis of auxin response using the R2D2 marker. Output is saved
in <path-to-repo>/processed_data.

"""

import os
project_path = '/home/henrik/projects/papers/bhatia_et_al_2019/'
os.chdir(os.path.join(project_path, 'code'))

import numpy as np
import mahotas as mh
import tifffile as tiff
from pycostanza.steepest import get_footprint, steepest_ascent
from pycostanza.labels import erode, dilate, merge_labels_distance, remove_labels_intensity
from pycostanza.labels import remove_labels_size, relabel, merge_labels_depth, merge_labels_small2closest
from scipy.ndimage.filters import gaussian_filter, median_filter
from skimage.morphology import remove_small_objects
from skimage.measure import regionprops
#from skimage.morphology import binary_closing, binary_opening

from misc import autocrop, listdir

### Import example inage
data_dir = os.path.join(project_path, 'intermediate_data', 'cropped_and_deconvolved_leaves')

files = listdir(data_dir, include='.tif')
files.sort()

def get_resolution(fobj):
    tags = fobj.pages.pages[0].tags
    x = float(tags['XResolution'].value[1]) / tags['XResolution'].value[0]
    y = float(tags['YResolution'].value[1]) / tags['YResolution'].value[0]
    z = fobj.imagej_metadata['spacing']
    return np.array([z, y, x])

for fname in files:
    ### Read in data
    sample, leaf, side = os.path.basename(fname).split('_')[0:3]
    sample = int(sample)
    f = tiff.TiffFile(fname)
    data = f.asarray()
    data = autocrop(data, fct=np.max, threshold=0)
    data = data.astype(np.float64)
    resolution = get_resolution(f)
    del f
    
    ### Use both channels for the segmentation
    int_img = data.copy()
    int_img[:, 1] /= np.max(int_img[:, 1])
    int_img[:, 2] /= np.max(int_img [:, 2])
    int_img = int_img [:, 1:]
    int_img = np.max(int_img, axis=1)
    int_img *= np.iinfo(np.uint16).max
    int_img = int_img.astype(np.uint16)

    # Preprocess by creating an initial binary mask
    footprint = get_footprint(3, 2)
    mask = int_img < 1. * mh.otsu(int_img, True)
    mask = remove_small_objects(mask, min_size=200, connectivity=2)
    #mask = binary_closing(mask, footprint)
    #mask = binary_opening(mask, footprint)
    int_img[mask] = 0

    ### Smooth the signal
    smooth_img = int_img.copy()
    smooth_img = median_filter(smooth_img, footprint=footprint)

    smooth_img = gaussian_filter(smooth_img, sigma=[.5, 1, 1])
    smooth_img = gaussian_filter(smooth_img, sigma=[.5, 1, 1])
    smooth_img = gaussian_filter(smooth_img, sigma=[.5*2/3, 1*2/3, 1*2/3])
    smooth_img = gaussian_filter(smooth_img, sigma=[.5*1/3, 1*1/3, 1*1/3])

    ### Perform initial segmentation
    lab_img = steepest_ascent(smooth_img, resolution, connectivity=2, mask=np.logical_not(mask))

    ### Merge and remove labels
    lab_img = merge_labels_distance(lab_img, smooth_img, threshold=1.5, resolution=resolution)
    lab_img = merge_labels_small2closest(lab_img, threshold=25, distance_upper_bound=3., 
                                         resolution=None) # 
    lab_img = remove_labels_size(lab_img, min_size=200, max_size=None, resolution=None)
    #        lab_img = merge_labels_depth(lab_img, int_img, threshold=2000., connectivity=2)
    lab_img = remove_labels_intensity(lab_img, int_img, threshold=3000.)
    lab_img = remove_labels_intensity(lab_img, data[:, 1], threshold=20.)

    lab_img = relabel(lab_img)

    ''' DATA COLLECTION '''
    ### Normalise to mDII frame-of-reference
    diimax = np.max(mh.labeled.labeled_max(data[:, 1], lab_img))
    mdiimax = np.max(mh.labeled.labeled_max(data[:, 2], lab_img))

    data[:, 1] *= diimax / mdiimax / mdiimax
    data[:, 2] /= mdiimax

    # Get expression / ratios
    rp1 = regionprops(lab_img, data[:, 1])
    rp2 = regionprops(lab_img, data[:, 2])
    ratios = [rp1[ii].mean_intensity / rp2[ii].mean_intensity for ii in xrange(len(rp1))]

    # Save to file    
    output_dir = os.path.join(project_path, 'processed_data', 'segmentation_data')
    output_path = os.path.join(output_dir, os.path.splitext(os.path.basename(fname))[0] + 
                               '_segmentation_data.dat')
    with open(output_path, 'w') as f:
        f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                'sample', 'leaf', 'side', 'label',
                'DII_mean', 'DII_median',  'DII_min',  'DII_max',
                'mDII_mean', 'mDII_median', 'mDII_min', 'mDII_max',
                'mean_ratio', 'median_ratio', 'min_ratio', 'max_ratio',
                'size',# 'minor_axis_length', 'major_axis_length',
                'solidity'))
        for ii in xrange(len(rp1)):
            f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    sample, leaf, side, rp1[ii].label,
                    
                    # DII
                    rp1[ii].mean_intensity, 
                    np.median(data[
                            rp1[ii].coords[:, 0], 
                            1,
                            rp1[ii].coords[:, 1], 
                            rp1[ii].coords[:, 2]]), 
                    rp1[ii].min_intensity, 
                    rp1[ii].max_intensity,
                    
                    # mDII
                    rp2[ii].mean_intensity, 
                        np.median(data[
                            rp2[ii].coords[:, 0], 
                            2,
                            rp2[ii].coords[:, 1], 
                            rp2[ii].coords[:, 2]]), 
                    rp2[ii].min_intensity, 
                    rp2[ii].max_intensity,

                    # Ratio
                    rp1[ii].mean_intensity / rp2[ii].mean_intensity, 
                    np.median(data[
                            rp1[ii].coords[:, 0],
                            1,
                            rp1[ii].coords[:, 1], 
                            rp1[ii].coords[:, 2]]) / (np.median(data[
                            rp2[ii].coords[:, 0], 
                            2,
                            rp2[ii].coords[:, 1], 
                            rp2[ii].coords[:, 2]]) + 1e-15), 
                    rp1[ii].min_intensity / (rp2[ii].min_intensity + 1e-15), 
                    rp1[ii].max_intensity / rp2[ii].max_intensity,
                    
                    # Structural
                    len(rp1[ii].coords) * np.product(resolution), # rp1[ii].minor_axis_length, rp1[ii].major_axis_length,
                    rp1[ii].solidity))

    # Save
    tiff.imsave(os.path.join(project_path, 'processed_data', 'segmented_leaves', 
                         os.path.splitext(os.path.basename(fname))[0] + '_segmented.tif'), 
        data=lab_img.astype(np.uint16), shape=lab_img.shape, dtype=np.uint16, imagej=True, 
        metadata={'spacing' : resolution[0]}, resolution=1. / resolution[1:])
