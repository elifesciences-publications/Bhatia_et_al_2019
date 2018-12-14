# Hello
This repository is made available for the publication "Quantitative analysis of auxin sensing in leaf primordia argues against proposed role in regulating leaf dorsoventrality",
submitted to *eLife* 18-06-2018, and authored by:

- Neha Bhatia (bhatia@mpipz.mpg.de)
- Henrik Åhl (henrik.aahl@slcu.cam.ac.uk)
- Henrik Jönsson (henrik.jonsson@slcu.cam.ac.uk)
- Marcus Heisler (marcus.heisler@sydney.edu.au)

For queries relating to the paper, contact Marcus Heisler (marcus.heisler@sydney.edu.au).
Questions related to the code are best addressed to Henrik Åhl (henrik.aahl@slcu.cam.ac.uk).

# Installing prerequisites
## Using pip
In order to run the code and replicate the output, install required packages by running:
```
pip install numpy pandas scipy tifffile scikit-image==0.14 pycostanza mahotas==1.4.5 matplotlib==2.2.3 scikit-learn==0.20.0
```
## Using Anaconda
If you are running Anaconda, you can find an environment file in [code/environment.yml](https://gitlab.com/slcu/teamHJ/publications/bhatia_et_al_2019/blob/master/code/environment.yml). Install the environment by running
```
conda env create -f code/environment.yml
```
from the repository root directory.

# Replicating data
The data can be replicated by running the scripts
```python
code/deconvolution.py
code/segmentation.py
code/figure_generation.py
```
in succession. Note that the deconvolution uses large amounts of RAM, which might prove troublesome for small-scale desktop computers.