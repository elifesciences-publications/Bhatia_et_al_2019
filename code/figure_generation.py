#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 14:15:20 2018

@author: henrik
"""

import os
project_path = '/home/henrik/projects/r2d2_leaves/'
os.chdir(os.path.join(project_path, 'code'))

import scipy
import numpy as np
import pandas as pd
import seaborn as sns
import tifffile as tiff
from scipy.stats import sem
from sklearn.metrics import auc
import matplotlib.pyplot as plt
from scipy.stats import normaltest
from pycostanza.labels import dilate
from skimage.measure import regionprops

from misc import listdir, get_resolution, natural_sort, match_shape

# Get paths
data_dir = os.path.join(project_path, 'processed_data', 'segmentation_data')
data_files = listdir(data_dir, include='.dat')
data_files.sort()
img_dir = os.path.join(project_path, 'processed_data', 'segmented_leaves')
img_files = listdir(img_dir, include='.tif')
img_files.sort()

# Read in all data to one frame
frame = pd.concat([pd.read_csv(file_, index_col=None, header='infer', sep='\t') for file_ in data_files])

# Filter
frame = frame[frame['size'] > 20]
frame = frame[frame['size'] < 200]
frame = frame[frame['mean_ratio'] <= 1.0]
#frame = frame[frame['mean_ratio'] > .05]
#frame = frame[frame['DII_max'] > 1000]

# Do some statistics that might be of interest
ab = frame[frame['side'] == 'ab']
ad = frame[frame['side'] == 'ad']
ad = ad.sort_values(by='mean_ratio', ascending=True)
ab = ab.sort_values(by='mean_ratio', ascending=True)

print(ab.shape[0])
print(ad.shape[0])


print normaltest(ab.mean_ratio).pvalue
print normaltest(ad.mean_ratio).pvalue
print np.mean(ab.mean_ratio), np.std(ab.mean_ratio), sem(ab.mean_ratio)
print np.mean(ad.mean_ratio), np.std(ad.mean_ratio), sem(ad.mean_ratio)
print scipy.stats.ks_2samp(ad.mean_ratio, ab.mean_ratio)
print scipy.stats.mannwhitneyu(ad.mean_ratio, ab.mean_ratio, use_continuity=True, alternative='two-sided')

# Plot things
###############################################################################
### Statistics etc.
###############################################################################
''' Set Seaborn plot style '''
sns.set(style='whitegrid', font_scale=2,
        rc={'axes.grid': False,
              'figure.figsize': (7., 5.)})
plt.style.use('seaborn-whitegrid')
palette = plt.get_cmap('Set1')

''' Generate ranked cell comparison plot '''
import matplotlib.patches as mpatches
patchList = []
legend_dict = { 'Adaxial' : palette(1), 'Abaxial' : palette(2)}
for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

plt.plot((ad.mean_ratio.values), color=palette(1), linewidth=4)
plt.plot((ab.mean_ratio.values), color=palette(2), linewidth=4)
plt.legend(handles=patchList, loc=2, ncol=1)
plt.xlabel('Ranked cell index')
plt.ylabel('Mean intensity ratio')

''' Generate AUC curve to compare distributions '''
ranked_data = frame.sort_values(by='mean_ratio')[['side', 'mean_ratio']]
side = ranked_data['side'].values.astype('str')
vals = ranked_data['mean_ratio'].values
a = np.zeros(len(side) + 1)
b = np.zeros(len(side) + 1)
for ii, s in enumerate(side, start=1):
    if s == 'ad':
        a[ii] = a[ii - 1] + 1
        b[ii] = b[ii - 1]
    else:
        a[ii] = a[ii - 1]
        b[ii] = b[ii - 1] + 1

a = a / max(a)
b = b / max(b)

sns.set(style='whitegrid', font_scale=2,
        rc={'axes.grid': False,
              'figure.figsize': (7., 5.)})
plt.style.use('seaborn-whitegrid')
palette = plt.get_cmap('Set1')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.step(b, a, color=palette(1), linewidth=2)
ax.plot(np.linspace(0, 1, 100), np.linspace(0, 1, 100), color='gray', linewidth=2, linestyle='--')
plt.xlabel('Abaxial ranked fraction')
plt.ylabel('Adaxial ranked fraction')
ax.text(.05,.9, 'AUC = {0:.3g}'.format(auc(b, a, True)), fontsize=24)
print auc(b, a, True)
ax.locator_params(axis='y', nbins=5)
ax.locator_params(axis='x', nbins=5)
fig.subplots_adjust(left=.15, bottom=.18)
fig.savefig(os.path.join(project_path, 'figures', 'AUC.png'))
fig.savefig(os.path.join(project_path, 'figures', 'AUC.pdf'))

''' Generate adaxial-abaxial overall distribution comparison plot -- MEAN '''
sns.set(style='whitegrid', font_scale=2,
        rc={'axes.grid': False,
              'figure.figsize': (7., 5.)})

plt.style.use('seaborn-whitegrid')
palette = plt.get_cmap('Set1')
ax = sns.violinplot(x="side", y="mean_ratio", data=frame,
                    inner='box', width=.8, linewidth=2, palette='muted', bw=.2, cut=0)
ax.set(xticklabels=['Abaxial', 'Adaxial'])
plt.xlabel('')
plt.ylabel('Mean intensity ratio')
plt.savefig(os.path.join(project_path, 'figures', 'adab_violins_mean.png'))
plt.savefig(os.path.join(project_path, 'figures', 'adab_violins_mean.pdf'))

''' Generate adaxial-abaxial overall distribution comparison plot -- MEDIAN '''
sns.set(style='whitegrid', font_scale=2,
        rc={'axes.grid': False,
              'figure.figsize': (7., 5.)})

plt.style.use('seaborn-whitegrid')
palette = plt.get_cmap('Set1')
ax = sns.violinplot(x="side", y="median_ratio", data=frame,
                    inner='box', width=.8, linewidth=2, palette='muted', bw=.2, cut=0)
ax.set(xticklabels=['Abaxial', 'Adaxial'])
plt.xlabel('')
plt.ylabel('Median intensity ratio')
plt.savefig(os.path.join(project_path, 'figures', 'adab_violins_median.png'))
plt.savefig(os.path.join(project_path, 'figures', 'adab_violins_median.pdf'))

''' Generate adaxial-abaxial overall distribution comparison plot -- MAX '''
sns.set(style='whitegrid', font_scale=2,
        rc={'axes.grid': False,
              'figure.figsize': (7., 5.)})

plt.style.use('seaborn-whitegrid')
palette = plt.get_cmap('Set1')
ax = sns.violinplot(x="side", y="max_ratio", data=frame,
                    inner='box', width=.8, linewidth=2, palette='muted', bw=.2, cut=0)
ax.set(xticklabels=['Abaxial', 'Adaxial'])
plt.xlabel('')
plt.ylabel('Maximum intensity ratio')
plt.savefig(os.path.join(project_path, 'figures', 'adab_violins_max.png'))
plt.savefig(os.path.join(project_path, 'figures', 'adab_violins_max.pdf'))

''' Generate adaxial-abaxial overall distribution comparison plot -- MIN '''
sns.set(style='whitegrid', font_scale=2,
        rc={'axes.grid': False,
              'figure.figsize': (7., 5.)})

plt.style.use('seaborn-whitegrid')
palette = plt.get_cmap('Set1')
ax = sns.violinplot(x="side", y="min_ratio", data=frame,
                    inner='box', width=.8, linewidth=2, palette='muted', bw=.2, cut=0)
ax.set(xticklabels=['Abaxial', 'Adaxial'])
plt.xlabel('')
plt.ylabel('Minimum intensity ratio')
plt.savefig(os.path.join(project_path, 'figures', 'adab_violins_min.png'))
plt.savefig(os.path.join(project_path, 'figures', 'adab_violins_min.pdf'))


''' Generate sample-sample distribution comparison figure '''
sns.set(style='whitegrid', font_scale=2,
        rc={'axes.grid': False,
              'figure.figsize': (12., 6.)})

plt.style.use('seaborn-whitegrid')
palette = plt.get_cmap('Set1')

ax = sns.violinplot(x="side", y="mean_ratio", hue='sample', data=frame,
                    inner=None, opacity=.5, cut=0, palette='muted', width=.8)
ax = sns.swarmplot(x='side', y='mean_ratio', hue='sample',
                   data=frame, size=1, linewidth=2, dodge=True)
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[0:10], labels[0:10], bbox_to_anchor=(1, 1),
           loc=2, borderaxespad=.1, title='Sample')
ax.set(xticklabels=['Abaxial', 'Adaxial'])
plt.xlabel('')
plt.ylabel('Mean intensity ratio')
plt.subplots_adjust(left=.1, bottom=.05, right=.85)
plt.savefig(os.path.join(project_path, 'figures', 'adab_violins_all.png'))
plt.savefig(os.path.join(project_path, 'figures', 'adab_violins_all.pdf'))

###############################################################################
### Spatial distributions
###############################################################################
for ii, fname in enumerate(img_files):
    f = tiff.TiffFile(fname)
    sample, leaf, side = os.path.basename(fname).split('_')[0:3]
    sample = int(sample)
    data = frame.loc[(frame['sample'] == sample) & (frame['leaf'] == leaf) & (frame['side'] == side)]
    ratios = (data.DII_mean / data.mDII_mean).values
    labels = data['label'].values
    
    lab_img = f.asarray()
    resolution = get_resolution(f)
    rp = regionprops(lab_img)
    rp = [rr for rr in rp if rr.label in labels]
    del f

    # Plot centroids
    ratio_plot = np.zeros(lab_img.shape)
    centroids = np.array([np.array(rr.centroid) for rr in rp])
    centroids = centroids.astype(np.uint16)
    ratio_plot[centroids[:,0], centroids[:,1], centroids[:,2]] = ratios
    
    tiff.imsave(os.path.join(project_path, 'figures', 'spatial', os.path.basename(fname)[:-14] + '_ratio_centroids.tif'), 
                dilate(np.max(ratio_plot, axis=0), size=5).astype(np.float32), dtype=np.float32, imagej=True,
                resolution=1. / resolution[1:])
    
    # Plot full nuclei
    ratio_plot = np.zeros(lab_img.shape)
    for jj, label in enumerate(data.label):
        ratio_plot[lab_img == label] = ratios[jj]
    
    tiff.imsave(os.path.join(project_path, 'figures', 'spatial', os.path.basename(fname)[:-14] + '_ratio_full.tif'), 
                ratio_plot.astype(np.float32), resolution=1. / resolution[1:], imagej=True, 
                metadata={'spacing' : resolution[0]})  
    
files = listdir('/home/henrik/projects/r2d2_leaves/figures/spatial', include='ratio_centroids')
files = natural_sort(files)
ab = files[::2]
ad = files[1::2]
add = [tiff.imread(ad[ii]) for ii in xrange(len(ab))]
abb = [tiff.imread(ab[ii]) for ii in xrange(len(ab))]

s1 = np.array([aa.shape for aa in add])
s2 = np.array([bb.shape for bb in abb])
max_shape = (np.max([s1[:, 0], s2[:, 0]]), np.max([s1[:, 1], s2[:, 1]]))

add = [match_shape(ss, max_shape) for ss in add]
abb = [match_shape(ss, max_shape) for ss in abb]

a = np.zeros((max_shape[0]*len(ab), max_shape[1]))
for ii in xrange(len(ab)):
    a[max_shape[0]*ii:(max_shape[0]*ii + abb[ii].shape[0]), :abb[ii].shape[1]] = abb[ii]

b = np.zeros((max_shape[0]*len(ab), max_shape[1]))
for ii in xrange(len(ab)):
    b[max_shape[0]*ii:(max_shape[0]*ii + add[ii].shape[0]), :add[ii].shape[1]] = add[ii]

#a = a.reshape(a.shape[0] / 2, a.shape[1] * 2)
#b = b.reshape(b.shape[0] / 2, b.shape[1] * 2)
    
samps = [aa.split('/spatial/')[1].split('_')[0] for aa in ad]
sides = ['L', 'R']*(len(ad)/2)
xlabs = [samps[ii] + sides[ii] for ii in xrange(len(samps))]

a1 = a[:(a.shape[0]/2)]
a2 = a[(a.shape[0]/2):]
b1 = b[:(b.shape[0]/2)]
b2 = b[(b.shape[0]/2):]

sns.set(style='whitegrid', font_scale=2,
        rc={'axes.grid': False,
              'figure.figsize': (8., 13.)})
fig = plt.figure()
bottom1 = 0.11
ax1 = fig.add_axes([0.1, bottom1, 0.2, 0.85])
ax2 = fig.add_axes([0.275, bottom1 , 0.2, 0.85])
ax3 = fig.add_axes([0.5, bottom1 , 0.2, 0.85])
ax4 = fig.add_axes([0.675, bottom1 , 0.2, 0.85])
cbar_ax = fig.add_axes([0.15, 0.08, 0.675, 0.01])

plt.style.use('seaborn-whitegrid')
palette = plt.get_cmap('viridis')
ax1 = sns.heatmap(a1, ax=ax1, vmax=1, cbar=False, xticklabels=[], yticklabels=[], mask=a1==0, cmap=palette)
ax2 = sns.heatmap(b1, ax=ax2, vmax=1, cbar=False, xticklabels=[], yticklabels=[], mask=b1==0, cmap=palette)
ax3 = sns.heatmap(a2, ax=ax3, vmax=1, cbar=False, xticklabels=[], yticklabels=[], mask=a2==0, cmap=palette)
ax4 = sns.heatmap(b2, ax=ax4, vmax=1, cbar=False, xticklabels=[], yticklabels=[], mask=b2==0, cmap=palette)

ax1.set_xlabel('Ad')    
ax1.xaxis.set_label_position('top')
ax2.set_xlabel('Ab')    
ax2.xaxis.set_label_position('top')
ax3.set_xlabel('Ad')    
ax3.xaxis.set_label_position('top')
ax4.set_xlabel('Ab')    
ax4.xaxis.set_label_position('top')

ax1.set_yticks(np.arange(b.T.shape[1] / len(ad) / 2, b.T.shape[1], b.T.shape[1] / len(ad))[:10])
ax1.set_yticklabels(xlabs[:(len(xlabs)/2)])
ax4.yaxis.tick_right()
ax4.yaxis.set_ticks_position('none') 
ax4.set_yticks(np.arange(b.T.shape[1] / len(ad) / 2, b.T.shape[1], b.T.shape[1] / len(ad))[:10])
ax4.set_yticklabels(xlabs[(len(xlabs)/2):])

import matplotlib as mpl
cmap = mpl.cm.cool
norm = mpl.colors.Normalize(vmin=0, vmax=1)
cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=palette, norm=norm,
                                orientation='horizontal')
#cb1.ax.set_xticks([0, 1], ['0', '1'])
cb1.ax.set_xticklabels([0, .2, .4, .6, .8, 1])
cb1.ax.set_label('')
#cb1.ax.set_label('Mean intensity ratio')
fig.text(.325, 0.03, 'Mean intensity ratio')
#cb1.set_ticks([])

fig.savefig(os.path.join(project_path, 'figures', 'spatial_distribution_centroids.png'))
fig.savefig(os.path.join(project_path, 'figures', 'spatial_distribution_centroids.pdf'))
