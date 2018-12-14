#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 17:43:39 2018

@author: henrik
"""

import numpy as np
import os
import re

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def get_resolution(fobj):
    tags = fobj.pages.pages[0].tags
    x = float(tags['XResolution'].value[1]) / tags['XResolution'].value[0]
    y = float(tags['YResolution'].value[1]) / tags['YResolution'].value[0]
    z = fobj.imagej_metadata['spacing']
    return np.array([z, y, x])

def autocrop(arr, threshold=8e7, fct=np.sum, channel=-1):
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

    return arr[cp[0, 0]:cp[0, 1], :, cp[1, 0]:cp[1, 1], cp[2, 0]:cp[2, 1]]


def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def listdir(path, include=None, exclude=None, full=True):
    files = os.listdir(path)
    files = np.array(files)

    if full:
        files = np.array(map(lambda x: os.path.join(path, x), files))

    # Include
    if isinstance(include, str):
        files = np.array(filter(lambda x: include in x, files))
    elif isinstance(include, (list, np.ndarray)):
        matches = np.array([np.array([inc in ii for ii in files]) for inc in include])
        matches = np.any(matches, axis=0)
        files = files[matches]

    # Exclude
    if isinstance(exclude, str):
        files = np.array(filter(lambda x: exclude not in x, files))
    elif isinstance(exclude, (list, np.ndarray)):
        matches = np.array([np.array([exc in ii for ii in files]) for exc in exclude])
        matches = np.logical_not(np.any(matches, axis=0))
        files = files[matches]

    return files

def match_shape(a, t, side='both', val=0):
    """

    Parameters
    ----------
    a : np.ndarray
    t : Dimensions to pad/trim to, must be a list or tuple
    side : One of 'both', 'before', and 'after'
    val : value to pad with
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
