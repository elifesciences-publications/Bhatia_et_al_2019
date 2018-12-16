# setup.py
# Usage: ``python setup.py build_ext --inplace``
from distutils.core import setup, Extension
import numpy
setup(name='_psf',
	ext_modules=[Extension('_psf', ['psf.c'],
	include_dirs=[numpy.get_include()])])
