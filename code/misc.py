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

This particular code file contains miscellaneous auxillary functions used by the
other provided scripts.

"""

import os
import re
import numpy as np

def natural_sort(l): 
    """ Perform natural sorting. """
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def get_resolution(fobj):
    """
    Get ZYX resolutions from a TiffFile object. Assuming the tags 'XResolution' and
    'YResolution' exist, as well as the metadata entry 'spacing'.
    
    """
 
    tags = fobj.pages.pages[0].tags
    x = float(tags['XResolution'].value[1]) / tags['XResolution'].value[0]
    y = float(tags['YResolution'].value[1]) / tags['YResolution'].value[0]
    z = fobj.imagej_metadata['spacing']
    return np.array([z, y, x])

def autocrop(arr, threshold=8e7, fct=np.sum, channel=-1):
    """
    Automatically crop a given (3D) image based on its intensity profile. 
    
    Parameters
    ----------
    arr : np.array
        Input image in ZCYX order.
    threshold : float
        Threshold setting where to crop to.
    fct : function 
        Function to use for operating over the slices.
    channel : int
        Which channel (index) to use for processing. If -1, the maximal value of all channels
        are used.
        
    Returns
    ----------
    arr : np.array
        The cropped array.
        
    
    """

    if channel == -1:
        sumarr = np.max(arr, axis=1)
    elif isinstance(channel, (list, np.ndarray, tuple)):
        sumarr = np.max(arr.take(channel, axis=1), axis=1)
    else:
        sumarr = sumarr[:, channel]

    cp = np.zeros((sumarr.ndim, 2), dtype=np.int)
    for ii in xrange(sumarr.ndim):
        summers = np.array([0, 1, 2])[np.array([0, 1, 2]) != ii]

        vals = fct(sumarr, axis=tuple(summers))
        first = next((e[0] for e in enumerate(vals) if e[1] > threshold), 0)
        last = len(
            vals) - next((e[0] for e in enumerate(vals[::-1]) if e[1] > threshold), 0)

        cp[ii] = first, last
    arr = arr[cp[0, 0]:cp[0, 1], :, cp[1, 0]:cp[1, 1], cp[2, 0]:cp[2, 1]]
    return arr 


def mkdir(path):
    """
    Create directory `path`, including necessary sub-directories, if not already 
    existing. Equivalent to the shell command mkdir -p <path>.

    """

    if not os.path.exists(path):
        os.makedirs(path)

def listdir(path, include=None, exclude=None, full=True):
    """
    List files in a given directory. 
    
    Parameters
    ----------
    path : str
        The directory path
    include : str or list of str
        Only include files mathing the given string(s).
    exclude : str or list of str
        Exclude files mathing the given string(s).
    full : bool
        If True, returns full path to file.
        
    Returns
    ----------
    files : np.array
        List of files in `path`.
    
    """
    
    files = os.listdir(path)
    files = np.array(files)

    if full:
        files = np.array(map(lambda x: os.path.join(path, x), files))

    # Include
    if isinstance(include, str):
        files = np.array(filter(lambda x: include in x, files))
    elif isinstance(include, (list, np.ndarray, tuple)):
        matches = np.array([np.array([inc in ii for ii in files]) for inc in include])
        matches = np.any(matches, axis=0)
        files = files[matches]

    # Exclude
    if isinstance(exclude, str):
        files = np.array(filter(lambda x: exclude not in x, files))
    elif isinstance(exclude, (list, np.ndarray, tuple)):
        matches = np.array([np.array([exc in ii for ii in files]) for exc in exclude])
        matches = np.logical_not(np.any(matches, axis=0))
        files = files[matches]

    return files

def match_shape(a, t, side='both', val=0):
    """
    Match an input array to a target shape by concatenating `val` to the given 
    array on the specified `side`. If the target array is smaller, the array is 
    trimmed using the same logic. 

    This function was a bit shamelessly stolen and modified from a (public) internet 
    source I can no longer find. If the author wants to make him- or herself known
    and recieve recognition, please contact henrik.aahl@slcu.cam.ac.uk.

    Parameters
    ----------
    a : np.ndarray
    t : Dimensions to pad/trim to, must be a list or tuple
    side : One of 'both', 'before', and 'after'
    val : value to pad with
    
    Returns
    ----------
    b : np.array
        Input array modified to fit the target shape of `t`.
    
    """
    try:
        if len(t) != a.ndim:
            raise TypeError(
                't shape must have the same number of dimensions as the input')
    except TypeError:
        raise TypeError('t must be array-like')

    try:
        if isinstance(val, (int, long, float, complex)):
            b = np.ones(t, a.dtype) * val
        elif val == 'max':
            b  = np.ones(t, a.dtype) * np.max(a)
        elif val == 'mean':
            b  = np.ones(t, a.dtype) * np.mean(a)
        elif val == 'median':
            b  = np.ones(t, a.dtype) * np.median(a)
        elif val == 'min':
            b  = np.ones(t, a.dtype) * np.min(a)
    except TypeError:
        raise TypeError('Pad value must be numeric or string')
    except ValueError:
        raise ValueError('Pad value must be scalar or valid string')

    aind = [slice(None, None)] * a.ndim
    bind = [slice(None, None)] * a.ndim

    # pad/trim comes after the array in each dimension
    if side == 'after':
        for dd in xrange(a.ndim):
            if a.shape[dd] > t[dd]:
                aind[dd] = slice(None, t[dd])
            elif a.shape[dd] < t[dd]:
                bind[dd] = slice(None, a.shape[dd])
    # pad/trim comes before the array in each dimension
    elif side == 'before':
        for dd in xrange(a.ndim):
            if a.shape[dd] > t[dd]:
                aind[dd] = slice(int(a.shape[dd] - t[dd]), None)
            elif a.shape[dd] < t[dd]:
                bind[dd] = slice(int(t[dd] - a.shape[dd]), None)
    # pad/trim both sides of the array in each dimension
    elif side == 'both':
        for dd in xrange(a.ndim):
            if a.shape[dd] > t[dd]:
                diff = (a.shape[dd] - t[dd]) / 2.
                aind[dd] = slice(int(np.floor(diff)), int(a.shape[dd] - np.ceil(diff)))
            elif a.shape[dd] < t[dd]:
                diff = (t[dd] - a.shape[dd]) / 2.
                bind[dd] = slice(int(np.floor(diff)), int(t[dd] - np.ceil(diff)))
    else:
        raise Exception('Invalid choice of pad type: %s' % side)

    b[bind] = a[aind]

    return b
