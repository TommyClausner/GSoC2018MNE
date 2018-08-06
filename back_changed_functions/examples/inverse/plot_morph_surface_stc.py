"""
=============================
Morph surface source estimate
=============================

This example demonstrates how to morph an individual subject **surface source
estimate** to a common reference space. It will be demonstrated using the
SourceMorph class. Pre-computed data will be morphed based on an affine
transformation on the surface vertices towards the reference vertices, which in
this example will be 'fsaverage'.

The transformation will be applied to the surface source estimate. The result
will be a plot showing the inflated surface representation of 'fsaverage',
overlaid with the morphed source estimate.
"""
# Author: Tommy Clausner <tommy.clausner@gmail.com>
#
# License: BSD (3-clause)
import os

from mne import read_source_estimate, SourceMorph
from mne.datasets import sample
from nilearn.plotting import plot_glass_brain

print(__doc__)

###############################################################################
# Setup paths
sample_dir_raw = sample.data_path()
sample_dir = os.path.join(sample_dir_raw, 'MEG', 'sample')
subjects_dir = os.path.join(sample_dir_raw, 'subjects')

fname_stc = os.path.join(sample_dir, 'sample_audvis-meg')

###############################################################################
# Load example data

# Read stc from file
stc = read_source_estimate(fname_stc, subject='sample')

###############################################################################
# Morph SourceEstimate

# Initialize SourceMorph for SourceEstimate
source_morph = SourceMorph(subject_from='sample',  # Default: None
                           subject_to='fsaverage',  # Default
                           subjects_dir=subjects_dir,  # Default: None
                           src=None,  # Default
                           spacing=5)  # Default

# Morph data
stc_fsaverage = source_morph(stc)

###############################################################################
# Plot results

surfer_kwargs = dict(
    hemi='lh', subjects_dir=subjects_dir,
    clim=dict(kind='value', lims=[8, 12, 15]), views='lateral',
    initial_time=0.09, time_unit='s', size=(800, 800),
    smoothing_steps=5)

brain = stc_fsaverage.plot(**surfer_kwargs)
brain.add_text(0.1, 0.9, 'Morphed to fsaverage', 'title', font_size=16)
