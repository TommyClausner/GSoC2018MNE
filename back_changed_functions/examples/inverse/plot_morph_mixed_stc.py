"""
===========================
Morph mixed source estimate
===========================
"""
import os.path as op

import mne
from mne.datasets import sample
from mne import setup_volume_source_space, setup_source_space
from mne import make_forward_solution
from mne.io import read_raw_fif
from mne.minimum_norm import make_inverse_operator, apply_inverse_epochs

# Set dir
data_path = sample.data_path()
subject = 'sample'
data_dir = op.join(data_path, 'MEG', subject)
subjects_dir = op.join(data_path, 'subjects')
bem_dir = op.join(subjects_dir, subject, 'bem')

# Set file names
fname_aseg = op.join(subjects_dir, subject, 'mri', 'aseg.mgz')

fname_model = op.join(bem_dir, '%s-5120-bem.fif' % subject)
fname_bem = op.join(bem_dir, '%s-5120-bem-sol.fif' % subject)

fname_raw = data_dir + '/sample_audvis_filt-0-40_raw.fif'
fname_trans = data_dir + '/sample_audvis_raw-trans.fif'
fname_cov = data_dir + '/ernoise-cov.fif'
fname_event = data_dir + '/sample_audvis_filt-0-40_raw-eve.fif'

# List of sub structures we are interested in. We select only the
# sub structures we want to include in the source space
labels_vol = ['Left-Amygdala',
              'Left-Thalamus-Proper',
              'Left-Cerebellum-Cortex',
              'Brain-Stem',
              'Right-Amygdala',
              'Right-Thalamus-Proper',
              'Right-Cerebellum-Cortex']

# Setup a surface-based source space
src = setup_source_space(subject, subjects_dir=subjects_dir,
                         spacing='oct6', add_dist=False)

# Setup a volume source space
# set pos=7.0 for speed issue
vol_src = setup_volume_source_space(subject, mri=fname_aseg,
                                    pos=7.0,
                                    bem=fname_model,
                                    volume_label=labels_vol,
                                    subjects_dir=subjects_dir)
# Generate the mixed source space
src += vol_src

# compute the fwd matrix
fwd = make_forward_solution(fname_raw, fname_trans, src, fname_bem,
                            mindist=5.0,  # ignore sources<=5mm from innerskull
                            meg=True, eeg=False,
                            n_jobs=1)

# Load data
raw = read_raw_fif(fname_raw, preload=True)
noise_cov = mne.read_cov(fname_cov)
events = mne.read_events(fname_event)

# Add a bad channel
raw.info['bads'] += ['MEG 2443']

# Pick MEG channels
picks = mne.pick_types(raw.info, meg=True, eeg=False, stim=False, eog=True,
                       exclude='bads')

# Define epochs for left-auditory condition
event_id, tmin, tmax = 1, -0.2, 0.5
epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
                    baseline=(None, 0), reject=dict(mag=4e-12, grad=4000e-13,
                                                    eog=150e-6))

# Compute inverse solution and for each epoch
snr = 1.0  # use smaller SNR for raw data
inv_method = 'dSPM'
parc = 'aparc'  # the parcellation to use, e.g., 'aparc' 'aparc.a2009s'

lambda2 = 1.0 / snr ** 2

# Compute inverse operator
inverse_operator = make_inverse_operator(raw.info, fwd, noise_cov,
                                         depth=None, fixed=False)

stcs = apply_inverse_epochs(epochs, inverse_operator, lambda2, inv_method,
                            pick_ori=None, return_generator=False)

labels_parc = mne.read_labels_from_annot(subject, parc=parc,
                                         subjects_dir=subjects_dir)
src = inverse_operator['src']

label_ts = mne.extract_label_time_course(stcs, labels_parc, src,
                                         mode='mean_flip',
                                         allow_empty=True,
                                         return_generator=False)

morph = mne.SourceMorph(subject_from=src[0]['subject_his_id'],
                        src=src, subjects_dir=subjects_dir)