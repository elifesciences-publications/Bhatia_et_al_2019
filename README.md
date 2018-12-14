# Hello
This repository contains material relating to tbe paper _Auxin is not asymmetrically distributed in initiating Arabidopsis leaves_ (Bhatia et al., 2018). 

# Installing prerequisites
## Using pip
In order to run the code and replicate the output, install required packages by running:
```
pip install numpy pandas scipy tifffile scikit-image==0.14 pycostanza mahotas==1.4.5 matplotlib==2.2.3 scikit-learn==0.20.0
```
## Using Anaconda
If you are running Anaconda, you can find an environment file in [code/environment.yml](https://gitlab.com/slcu/teamHJ/running_projects/r2d2_leaves/blob/master/code/environment.yml). Install the environment by running
```
conda env create -f code/environment.yml
```
from the repository root directory.

# Replicating data
The data can be replicated by running the scripts
```python
/code/deconvolution.py
/code/segmentation.py
/code/figure_generation.py
```
in succession. Note that the deconvolution uses large amounts of RAM, which might prove troublesome for small-scale desktop computers.