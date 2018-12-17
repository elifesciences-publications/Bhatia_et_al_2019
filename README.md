# Introduction
This repository is made available for the publication "Quantitative analysis of auxin sensing in leaf primordia argues against proposed role in regulating leaf dorsoventrality",
submitted to *eLife* 18-06-2018, and authored by:

- Neha Bhatia (bhatia@mpipz.mpg.de)
- Henrik Åhl (henrik.aahl@slcu.cam.ac.uk)
- Henrik Jönsson (henrik.jonsson@slcu.cam.ac.uk)
- Marcus Heisler (marcus.heisler@sydney.edu.au)

For queries relating to the paper, contact Marcus Heisler (marcus.heisler@sydney.edu.au).
Questions related to the code are best addressed to Henrik Åhl (henrik.aahl@slcu.cam.ac.uk).

# Installing prerequisites
## Custom packages
### Scikit-Image
We utilized Scikit-Image (0.14.1) for several of the image quantification steps. 
However, the distribution as-is contains a bug which causes problems when deconvolving 
images which contain values close or equal to 0. We therefore corrected this part 
of the code by using the machine epsilon for close-to-zero values. Specifically, 
the change consists of exchanging the code snippet 

```python
relative_blur = image / convolve_method(im_deconv, psf, 'same')
```

which is found on on line 389 in *skimage/restoration/deconvolution.py*, with

```python 
eps = np.finfo(image.dtype).eps
x = convolve_method(im_deconv, psf, 'same')
np.place(x, x==0, eps)
relative_blur = image / x + eps
```

A modified version ready for install is provided in *code/external/scikit-image-0.14.1*.
For installation:
```bash
cd code/external/scikit-image-0.14.1
python setup.py install
cd -
```

### PyCostanza
For installation:
```bash
cd code/external/pycostanza-0.1.3
python setup.py install
cd -
```

PyCostanza can also be installed via pip by running
```bash
pip install pycostanza==0.1.3
```

### PSF
For installation:
```bash
cd code/external/psf-2018.02.07
python setup.py build_ext --inplace
python setup.py install
cd -
```

## Additional packages
In order to run the code and replicate the output, install required packages by running:
```
pip install numpy pandas scipy tifffile mahotas==1.4.5 matplotlib==2.2.3 scikit-learn==0.20.0
```

If you are running Anaconda, you can find an environment file in [code/environment.yml](https://gitlab.com/slcu/teamHJ/publications/bhatia_et_al_2019/blob/master/code/environment.yml). 
Install the environment by running

```bash
conda env create -f code/environment.yml
```

from the repository root directory. The custom packages will have to be installed 
as described above. 

# Reproducing data and figures
The data can be replicated by running the scripts

```python
code/deconvolution.py
code/segmentation.py
code/figure_generation.py
```

in succession. Note that the deconvolution uses large amounts of RAM, which might 
prove troublesome for small-scale desktop computers.