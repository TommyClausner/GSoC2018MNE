# -*- coding: utf-8 -*-
r"""
===================================
Background information on morphing
===================================

Here we give some background information on morphing in general,
and how it is done in MNE-Python in particular. Recommended reading and
corresponding information can be found in Gramfort *et al.* 2013 [1]_ and
Avants *et al.* 2009 [2]_ as well as in this `dipy example`_. For a shorter and
more applied version of this tutorial see
:ref:`sphx_glr_auto_tutorials_plot_morph.py`.

.. contents::
    :local:

Problem statement
=================

Modern neuro imaging techniques such as souce reconstruction or fMRI analyses,
make use of advanced mathematical models and hardware to map brain activation
patterns into a subject specific anatomical brain space.
This enables the study of spatio-temporal brain activity. Amongst many others,
the representation of spatio-temporal activation patterns is often done by
overlaying the actual anatomical brain structure, with the respective
activation at the respective anatomical location. Hereby volumetric anatomical
MR images are often used as such or are transformed into an inflated surface.

It becomes obvious that in order to compute group level statistics, data
representations across subjects must be morphed to a common frame, such that
anatomically / functional similar structures are represented at the same
spatial location.

.. _tut_morphing_basics:

Morphing basics
================

Morphing describes the procedure of transforming a data representation in n-
dimensional space into another data representation in the same space. In the
context of neuroimaging data this space will mostly be 3-dimensional and
necessary to bring individual subject spaces into a common frame.
In general morphing operations can be split into two different kinds of
transformation: linear and non-linear morphs or mappings.

A mapping is linear if it satisfies the following two conditions:

.. math::    f(u + v) = f(u) + f(v)\ ,
.. math::    f(cu) = cf(u)\ ,

where :math:`u` and :math:`v` are from the same vector space and :math:`c` can
be any scalar. This means that any linear transform is a mixture of additive
and multiplicative operations and hence is often represented in terms of a
transformation matrix.

In turn, a non-linear mapping can include linear components, but furthermore
functions that are not limited by the above constraints. However it needs to be
understood that including non-linear operations will alter the relationship of
data points within a vector and thus cannot be represented as a transformation
matrix. Instead every data point can be mapped independent of other data
points. This becomes especially handy when morphing volumetric brain data. To
achieve a mapping of structurally similar areas between subjects, it is
inevitable to employ non-linear operations to account for morphological
differences that cannot be represented by a linear transform.

In MNE-Python "brain space" data is represented as source estimate data,
obtained by one of the implemented source reconstruction methods. See
:ref:`sphx_glr_auto_tutorials_plot_mne_dspm_source_localization.py`

It is thus represented as :class:`SourceEstimate <mne.SourceEstimate>`,
:class:`VectorSourceEstimate <mne.VectorSourceEstimate>`,
:class:`VolSourceEstimate <mne.VolSourceEstimate>` or a mixture of those.

The data in the first two cases is represented as "surfaces". This means that
the data is represented as vertices on an inflated brain surface
representation. In the last case the data is represented in a 4-dimensional
space, were the last dimension refers to the data's sample points.

Computing an inflated surface representation can be accomplished using
FreeSurfer. Thereby, spherical morphing of the surfaces can be employed to
bring data from different subjects into a common anatomical frame. When dealing
with volumetric data, a non-linear morph map - or better morph volume - based
on two different anatomical reference MRIs can be computed.

.. _tut_surface_morphing

Surface morphing
================

The spherical morphing of the surfaces accomplished by FreeSurfer can be
employed to bring data from different subjects into a common anatomical
frame. This chapter describes utilities which make use of the spherical
morphing procedure. mne_morph_labels morphs label files between subjects
allowing the definition of labels in a one brain and transforming them to
anatomically analogous labels in another. mne_average_estimates offers
the capability to compute averages of data computed with the MNE software
across subjects.

.. _tut_surface_morphing_maps:

Surface morphing - morph maps
=============================

The MNE software accomplishes surface morphing with help of morphing
maps which can be either computed on demand or precomputed using
mne_make_morph_maps ,see :ref:`tut_surface_morphing_precompute`. The morphing
is performed with help of the registered spherical surfaces (``lh.sphere.reg``
and ``rh.sphere.reg`` ) which must be produced in FreeSurfer .
A morphing map is a linear mapping from cortical surface values
in subject A (:math:`x^{(A)}`) to those in another
subject B (:math:`x^{(B)}`)

.. math::    x^{(B)} = M^{(AB)} x^{(A)}\ ,

where :math:`M^{(AB)}` is a sparse matrix
with at most three nonzero elements on each row. These elements
are determined as follows. First, using the aligned spherical surfaces,
for each vertex :math:`x_j^{(B)}`, find the triangle :math:`T_j^{(A)}` on the
spherical surface of subject A which contains the location :math:`x_j^{(B)}`.
Next, find the numbers of the vertices of this triangle and set
the corresponding elements on the *j* th row of :math:`M^{(AB)}` so that
:math:`x_j^{(B)}` will be a linear interpolation between the triangle vertex
values reflecting the location :math:`x_j^{(B)}` within the triangle
:math:`T_j^{(A)}`.

It follows from the above definition that in general

.. math::    M^{(AB)} \neq (M^{(BA)})^{-1}\ ,

*i.e.*,

.. math::    x_{(A)} \neq M^{(BA)} M^{(AB)} x^{(A)}\ ,

even if

.. math::    x^{(A)} \approx M^{(BA)} M^{(AB)} x^{(A)}\ ,

*i.e.*, the mapping is *almost* a
bijection.

.. _tut_surface_morphing_smoothing:

Surface morphing - smoothing
============================

The current estimates are normally defined only in a decimated
grid which is a sparse subset of the vertices in the triangular
tessellation of the cortical surface. Therefore, any sparse set
of values is distributed to neighboring vertices to make the visualized
results easily understandable. This procedure has been traditionally
called smoothing but a more appropriate name
might be smudging or blurring in
accordance with similar operations in image processing programs.

In MNE software terms, smoothing of the vertex data is an
iterative procedure, which produces a blurred image :math:`x^{(N)}` from
the original sparse image :math:`x^{(0)}` by applying
in each iteration step a sparse blurring matrix:

.. math::    x^{(p)} = S^{(p)} x^{(p - 1)}\ .

On each row :math:`j` of the matrix :math:`S^{(p)}` there
are :math:`N_j^{(p - 1)}` nonzero entries whose values
equal :math:`1/N_j^{(p - 1)}`. Here :math:`N_j^{(p - 1)}` is
the number of immediate neighbors of vertex :math:`j` which
had non-zero values at iteration step :math:`p - 1`.
Matrix :math:`S^{(p)}` thus assigns the average
of the non-zero neighbors as the new value for vertex :math:`j`.
One important feature of this procedure is that it tends to preserve
the amplitudes while blurring the surface image.

Once the indices non-zero vertices in :math:`x^{(0)}` and
the topology of the triangulation are fixed the matrices :math:`S^{(p)}` are
fixed and independent of the data. Therefore, it would be in principle
possible to construct a composite blurring matrix

.. math::    S^{(N)} = \prod_{p = 1}^N {S^{(p)}}\ .

However, it turns out to be computationally more effective
to do blurring with an iteration. The above formula for :math:`S^{(N)}` also
shows that the smudging (smoothing) operation is linear.

.. _tut_surface_morphing_precompute:

Surface morphing - precomputing
===============================

The utility mne_make_morph_maps was
created to assist mne_analyze and mne_make_movie in
morphing. Since the morphing maps described above take a while to
compute, it is beneficial to construct all necessary maps in advance
before using mne_make_movie .
The precomputed morphing maps are located in ``$SUBJECTS_DIR/morph-maps`` .
mne_make_morph_maps creates this directory automatically if it does not exist.
If this directory exists when mne_analyze or mne_make_movie is run
and morphing is requested, the software first looks for already
existing morphing maps there. Also, if mne_analyze or mne_make_movie have
to recompute any morphing maps, they will be saved to
``$SUBJECTS_DIR/morph-maps`` if this directory exists.

The names of the files in ``$SUBJECTS_DIR/morph-maps`` are
of the form:

 <*A*> - <*B*> -``morph.fif`` ,

where <*A*> and <*B*> are
names of subjects. These files contain the maps for both hemispheres,
and in both directions, *i.e.*, both :math:`M^{(AB)}` and :math:`M^{(BA)}`, as
defined above. Thus the files <*A*> - <*B*> -``morph.fif`` or <*B*> - <*A*> -
``morph.fif`` are functionally equivalent. The name of the file produced by
mne_analyze or mne_make_movie depends on the role of <*A*> and <*B*> in
the analysis.

If you choose to compute the morphing maps in batch in advance,
use :ref:`mne_make_morph_maps`.

.. _tut_volumetric_morphing:

Volumetric morphing
===================

The key difference in volumetric morphing as compared to morphing surfaces, is
that the data is represented as a volume. A volume is a 3-dimensional data
representation. It is necessary to understand that the difference to a mesh
(what's commonly meant when referring to "3D model") is that the mesh is
"empty", while the volume is not. Whereas the mesh is defined by the
vertices of the outer hull, the volume is defined by that and the points it is
containing. Hence morphing volumetric data does not only require to map the
surface data, but also it's content in the correct way.

It becomes more easy to understand, when thinking about morphing one brain to
another. We not only want the cortices to overlap anatomically as good as
possible, but also all sub-cortical structures.

In MNE-Python the implementation of volumetric morphing is achieved by wrapping
the corresponding functionality from DiPy. See this `dipy example`_ for
reference.

The volumetric morphed is implemented as a two stage approach. First two
reference brains are are aligned using an affine linear registration and later
a non-linear Symmetric Diffeomorphic Registration in 3D.

.. _tut_volumetric_morphing_affine:

Volumetric morphing - affine registration
=========================================

See `http://nipy.org/dipy/examples_built/affine_registration_3d.html` for
reference.

Our goal is to pre-align both reference volumes as good as possible, to make it
easier for the later non-linear optimization to converge to an acceptable
minimum. The quality of the pre-alignment will be assessed using the
mutual information that is aimed to be maximal [3]_.

Mutual information can be defined as the amount of predictive value two
variables share. Thus how much information about a random variable B can be
obtained through variable A.

It can further be expressed in terms of entropy as the difference between the
joined entropy of A and B the respective conditional entropies. Hence the
higher the joint entropy, the lower the conditional and hence one variable is
more predictive for the respective other. Aiming for maximizing the mutual
information can thus be seen as reducing the conditional entropy and thus the
amount of information required from a second variable, to describe the system.

If we find a transformation such, that both volumes are overlapping as good as
possible, then the location of a particular area in one brain, would be highly
predictive for the same location on the second brain. In turn mutual
information is high, whereas the conditional entropy is low.

The specific optimization algorithm used for mutual information driven affine
registration is described in Mattes *et al.* 2003 [4]_.

In essence, a gradient decent is used to minimize the negative mutual
information, while optimizing the set of parameters of the image discrepancy
function.

.. _tut_volumetric_morphing_sdr:

Volumetric morphing - Symmetric Diffeomorphic Registration
==========================================================

See `http://nipy.org/dipy/examples_built/syn_registration_3d.html` for
reference.

Symmetric Diffeomorphic Image Registration is described in Avants *et al.* 2009
[2]_.

A map between two objects (manifolds that need to be differentiable) is
diffeomorphic if it is invertable (so is it's inverse). Hence it can be seen
as a locally linear, smooth map, that describes how each point on one object
relates to the same point on a second object. Imagine a scrambled and an intact
sheet of paper. There is a clear mapping between each point of the first, to
each point of the second object.

The introduced "symmetry" refers to symmetry after implementation. That is that
morphing A to B yields computationally the same result as morphing B to A.

As optimization criterion the cross-correlation was chosen, which is a
description of similarity between two data series.

.. _tut_sourcemorph

:class:`SourceMorph <mne.SourceMorph>`
======================================

:class:`SourceMorph <mne.SourceMorph>` is MNE-Python's source estimation
morphing operator. It can perform all necessary computations, to achieve the
above transformations on surface source estimate representations as well as for
volumetric source estimates. This includes
:class:`SourceEstimate <mne.SourceEstimate>` and
:class:`VectorSourceEstimate <mne.VectorSourceEstimate>` for surface
representations and :class:`VolSourceEstimate <mne.VolSourceEstimate>` for
volumetric representations.

SourceMorph can take general a type specific keyword arguments. The following
general keyword arguments are accepted:

    * *subject_from*: string pointing to the respective subject folder
      representing the subject the is going to be morphed.
      E.g. subject_from='Bert'. The default is subject_from=None. Within the
      respective subject folder (e.g. SUBJECTS_DIR/Bert), the result of a
      anatomical segmentation, done with FreeSurfer should be present. More
      specifically FreeSurfer surface data is employed to achieve surface morph
      operations, whereas for volume source estimates brain.mgz will be used by
      default.
    * *subject_to*: string pointing to the respective subject folder
      representing the subject the is used as reference for the respective
      morph. Hence this it represents the target space. Similar data as for
      subject_from is used. The default is subject_to='fsaverage', hence is not
      otherwise specified fsaverage will be the target space. Note that it is
      possible as well to specify a path pointing to any brain.mgz that is used
      as reference volume.
    * *subjects_dir*: FreeSurfer subject directory. If not defined in
      SUBJECTS_DIR, subjects_dir should be define, otherwise the default is
      None, utilizing the path set in the environment.
    * *src*: The list of source space corresponding to the source estimate that
      is targeted as a morph. While for VolSourceEstimates, src is must be set,
      it is not necessary to be set when attempting to perform a surface morph.
      Hence the default is src=None. If no src is provided (only possible for
      surface morphs), then the SourceMorph operator is only set up and the
      computation of the morph takes place, once the object is called. In the
      case of volumetric data, src must be provided, in order to obtain the
      relevant information of how the data maps to anatomical data, that is
      used to compute the morph.
    * *spacing*: This parameter is fundamentally different depending on the
      underlying morph. Please see :ref:`tut_sourcemorph_surf` for information
      related to surface morphs and :ref:`tut_sourcemorph_vol` for information
      related to volumetric morphs. Default is spacing=5
    * *precomputed*: If precomputed morph data is already present, or custom
      made morph data will be used, it can be set by the keyword argument
      precomputed. If not None, the argument must be a dictionary, carrying
      type specific morphing information set up in the same way as used by
      SourceMorph. Please see :ref:`tut_sourcemorph_surf` and
      :ref:`tut_sourcemorph_vol` for information about type specific morph
      parameters. precomputed=morph_params will thus override my_morph.params
      The default is precomputed=None

A SourceMorph object can be created, by initializing an instance and setting
desired key word arguments like so: ``my_morph = mne.SourceMorph(...)``

``my_morph`` will have all arguments set or their default values as attributes.
Furthermore it indicates of which kind it is ``my_morph.kind`` and what are the
respective morphing parameters ``my_morph.params``.

``my_morph.params`` is a dictionary, that varies depending on the type of
morph, all relevant information is stored. See :ref:`tut_sourcemorph_surf` and
:ref:`tut_sourcemorph_vol` for information about type specific morph
parameters.

:ref:`tut_surface_morphing`.

.. _tut_sourcemorph_surf

:class:`SourceMorph <mne.SourceMorph>` for :class:`mne.SourceEstimate`
======================================================================

In addition to general keyword arguments, SourceMorph can take multiple
arguments depending on the underlying morph. For (Vector)SourceEstimates, those
keyword arguments include:

    * *spacing*: In case of (Vector)SourceEstimate spacing can be an integer a
      list of 2 np.array or None. The defaut is spacing=5. Spacing refers to
      what was known as grade in previous versions of MNE-Python. It defines
      the esolution of the icosahedral mesh (typically 5). If None, all
      vertices will be used (potentially filling the surface). If a list,
      then values will be morphed to the set of vertices specified in
      in spacing[0] and spacing[1]. Note that specifying the vertices (e.g.,
      grade=[np.arange(10242), np.arange(10242)] for fsaverage on a
      standard grade 5 source space) can be substantially faster than
      computing vertex locations. Note that if subject='fsaverage'
      and 'spacing=5', this set of vertices will automatically be used
      (instead of computed) for speed, since this is a common morph.
    * *smooth*: Number of iterations for the smoothing of the surface data.
      If None, smooth is automatically defined to fill the surface with
      non-zero values. The default is None.
    * *xhemi*: If True data can be morphed between hemispheres by setting. The
      full cross-hemisphere morph matrix maps left to right and right to left.
      A matrix for cross-mapping only one hemisphere can be constructed by
      specifying the appropriate vertices, for example, to map the right
      hemisphere to the left:
      ``vertices_from=[[], vert_rh], vertices_to=[vert_lh, []]``.
      Cross-hemisphere mapping requires appropriate ``sphere.left_right``
      morph-maps in the subject's directory. These morph maps are included
      with the ``fsaverage_sym`` FreeSurfer subject, and can be created for
      other subjects with the ``mris_left_right_register`` FreeSurfer command.
      The ``fsaverage_sym`` subject is included with FreeSurfer > 5.1 and can
      be obtained as described
      `here <http://surfer.nmr.mgh.harvard.edu/fswiki/Xhemi>`_. For statistical
      comparisons between hemispheres, use of the symmetric ``fsaverage_sym``
      model is recommended to minimize bias [5]_.

.. _tut_sourcemorph_vol

:class:`SourceMorph <mne.SourceMorph>` for :class:`mne.VolSourceEstimate`
=========================================================================

In addition to general keyword arguments, SourceMorph can take multiple
arguments depending on the underlying morph. For VolSourceEstimates, those
keyword arguments include:

    * *spacing*: In case of VolSourceEstimate spacing can be an integer, float,
      tuple of integer or float or None. The default is spacing=5. Spacing
      refers to the voxel size that is used to compute the volumetric morph.
      Since two volumes are compared "point wise" the number of slices in each
      orthogonal direction has direct influence on the computation time and
      accuracy of the morph. See :ref:`tut_volumetric_morphing_sdr` to
      understand why this is the case. Spacing thus can also be seen as the
      voxel size to which both reference volumes will be resliced before
      computing the symmetric diffeomorphic volume. An integer or float value,
      will be interpreted as isotropic voxel size in mm. Setting a tuple allows
      for anisotropic voxel sizes e.g. (1., 1., 1.2). If None the full
      resolution of the MRIs will be used. Note, that this can cause long
      computation times.
    * *niter_affine*: As described in :ref:`tut_volumetric_morphing_affine` an
      iterative process is used to find the transformation that maps one image
      to another. This iterative process is performed in multiple levels and a
      number of iterations per level. A level is a stage of iterative
      refinement with a certain level of precision. The higher or later the
      level the more refined the iterative optimization will be, requiring more
      computation time. The number of levels and the number of iterations per
      level are defined as a tuple of integers, where the number of integers or
      the length of the tuple defines the number of levels, whereas the integer
      values themselves represent the number of iterations in that respective
      level. The default is niter_affine=(100, 100, 10) referring to a 3 stage
      optimization using 100, 100 and 10 iterations for the 1st, 2nd and 3rd
      level.
    * *niter_sdr*: As described in :ref:`tut_volumetric_morphing_sdr` an
      iterative process is used to find the transformation that maps one image
      to another. This iterative process is performed in multiple levels
      similar to the affine optimization
      (:ref:`tut_volumetric_morphing_affine`). The default is
      niter_sdr=(5, 5, 3) referring to a 3 stage optimization using 5, 5 and 3
      iterations for the 1st, 2nd and 3rd level.

.. _tut_sourcemorph_methods

:class:`SourceMorph <mne.SourceMorph>` methods
==============================================

Once an instance of SourceMorph was created, it exposes 3 methods:

    * *:meth`my_morph() <mne.SourceMorph.__call__>`:*: Calling an instance of
      SourceMorph on :class:`SourceEstimate <mne.SourceEstimate>`,
      :class:`VectorSourceEstimate <mne.VectorSourceEstimate>` or
      :class:`VolSourceEstimate <mne.VolSourceEstimate>`, will apply the
      precomputed morph to the input data and return the morphed source
      estimate (``stc_morphed = my_morph(stc)``). If a surface morph was
      attempted and no :class:`source space <mne.SourceSpaces>` was provided
      during instantiation of SourceMorph, then the actual computation of the
      morph will take place, using the input data as reference data, rather
      then precomputing it based on the source space data.
      Additionally the method takes the same keyword arguments as
      :meth:`my_morph.as_volume() <mne.SourceMorph.as_volume>`, given that
      `as_volume=True`. This means that the result will not be a source
      estimate, but instead NIfTI image representing the source estimate data
      in the specified way. If `as_volume=False` all other volume related
      arguments will be ignored.
    * *:meth`my_morph.as_volume() <mne.SourceMorph.as_volume>`:*: This method
      only works with :class:`VolSourceEstimate <mne.VolSourceEstimate>`s. It
      returns a NIfTI image of the source estimate. *mri_resolution* can be
      defined to change the resolution of the output image.
      ``mri_resolution=True`` will output an image in the same resolution as
      the MRI information stored in :class:`src <mne.SourceSpaces>`. If
      ``mri_resolution=False`` the output image will have the same resolution
      as defined in 'spacing' when instantiating the morph (see
      :ref:`tut_sourcemorph_vol`). Furthermore, mri_resolution can be defined
      as integer, float or tuple of integer or float to refer to the desired
      voxel size in mm. A single value will be interpreted as isotropic voxel
      size, whereas anisotropic dimensions can be defined using a tuple. Note,
      that a tuple must have a length of 3 referring to the 3 orthogonal
      spatial dimensions. The default is mri_resolution=False.
      The keyword argument *mri_space* asks, whether to use the same reference
      space as the reference MRI of the reference space of the source estimate.
      The default is mri_space=True.
      Furthermore a keyword argument called *apply_morph* can be set,
      indicating whether to apply the precomputed morph. In combination with
      the keyword argument 'as_volume', this can be used to produce morphed and
      unmorphed NIfTIs. The default is apply_morph=False.
    * *:meth`my_morph.save() <mne.SourceMorph.save>`:*: Saves the morph object
      to disk. The only input argument is the filename. Since the object is
      stored in HDF5 ('.h5') format, the filename will be extended by
      '-morph.h5' if no file extension was initially defined. To read saved
      SourceMorph objects, use :func:`mne.read_source_morph`.

.. _tut_sourcemorph_alternative

Alternative API
===============

Some operations can be performed using the respective source estimate itself.
This is mostly to support the API of previous versions of MNE-Python.

In this tutorial we will morph different kinds of source estimation results
between individual subject spaces using :class:`mne.SourceMorph`.
For group level statistical analyses subject specific results have to be mapped
to a common space.

We will use precomputed data and morph surface and volume source estimates to a
common space. The common space of choice will be FreeSurfer's "fsaverage".

Furthermore we will convert our volume source estimate into a NIfTI image using
:meth:`morph.as_volume <mne.SourceMorph.as_volume>`.
"""
# Author: Tommy Clausner <tommy.clausner@gmail.com>
#
# License: BSD (3-clause)

import os

###############################################################################
# Setup
# -----
#
# We first import the required packages and define a list of filenames for
# various datasets we are going to use to run this tutorial.
import matplotlib.pylab as plt
import nibabel as nib
from mne import (read_evokeds, SourceMorph, read_source_estimate)
from mne.datasets import sample
from mne.minimum_norm import apply_inverse, read_inverse_operator
from nilearn.image import index_img
from nilearn.plotting import plot_glass_brain

###############################################################################
# Background
# ----------
#
# - why -> averaging
# - how surface (https://www.frontiersin.org/articles/10.3389/fnins.2013.00267/full)
# - how volume (difficulty with volumes, morphing brain T1s, get affine + sdr,
# don't forget optimization procedures)
# - more info somewhere else (sphinx gallery)

# We use the MEG and MRI setup from the MNE-sample dataset
sample_dir_raw = sample.data_path()
sample_dir = sample_dir_raw + '/MEG/sample'
subjects_dir = sample_dir_raw + '/subjects'

fname_evoked = sample_dir + '/sample_audvis-ave.fif'

fname_surf = os.path.join(sample_dir, 'sample_audvis-meg')
fname_vol = os.path.join(sample_dir,
                         'sample_audvis-grad-vol-7-fwd-sensmap-vol.w')

fname_inv_surf = os.path.join(sample_dir,
                              'sample_audvis-meg-eeg-oct-6-meg-eeg-inv.fif')
fname_inv_vol = os.path.join(sample_dir,
                             'sample_audvis-meg-vol-7-meg-inv.fif')

fname_t1_fsaverage = subjects_dir + '/fsaverage/mri/brain.mgz'

###############################################################################
# Data preparation
# ----------------
#
# First we load the respective example data for surface and volume source
# estimates
stc_surf = read_source_estimate(fname_surf, subject='sample')

# Afterwards we load the corresponding source spaces as well
src_surf = read_inverse_operator(fname_inv_surf)['src']
inv_src = read_inverse_operator(fname_inv_vol)
src_vol = inv_src['src']

# ensure subject is not None
src_vol[0]['subject_his_id'] = 'sample'

# For faster computation we redefine tmin and tmax
stc_surf.crop(0.09, 0.1)

evoked = read_evokeds(fname_evoked, condition=0, baseline=(None, 0))

# Apply inverse operator
stc_vol = apply_inverse(evoked, inv_src, 1.0 / 3.0 ** 2, "dSPM")

# To save memory
stc_vol.crop(0.09, 0.1)

###############################################################################
# Setting up :class:`mne.SourceMorph` for :class:`mne.SourceEstimate`
# -------------------------------------------------------------------
#
# - point out surface specific stuff from above and relate to code


# SourceMorph initialization If src is not provided, the morph will not be
# pre-computed but instead will be prepared for morphing when calling. This
# works only with (Vector)SourceEstimate

morph_surf = SourceMorph(subject_from='sample',  # Default: None
                         subject_to='fsaverage',  # Default
                         subjects_dir=subjects_dir,  # Default: None
                         src=None,  # Default
                         spacing=5,  # Default
                         smooth=None,  # Default
                         xhemi=False)  # Default

###############################################################################
# Setting up :class:`mne.SourceMorph` for :class:`mne.VolSourceEstimate`
# ----------------------------------------------------------------------
#
# - point out volume specific stuff from above and relate to code

# Ideally subject_from can be inferred from src, subject_to is 'fsaverage' by
# default and subjects_dir is set in the environment. In that case SourceMorph
# can be initialized taking only src as parameter.

morph_vol = SourceMorph(subject_from='sample',  # Default: None
                        subject_to='fsaverage',  # Default
                        subjects_dir=subjects_dir,  # Default: None
                        spacing=(3., 3., 3.),  # grid spacing (3., 3., 3.) mm
                        src=src_vol,  # Default: None
                        niter_affine=(100, 100, 10),  # Default
                        niter_sdr=(5, 5, 3))  # Default

###############################################################################
# In a nutshell
# -------------
#
# For many applications the respective morph will probably look like this:

# Compute morph
# morph = SourceMorph(subject_from='sample',  # Default: None
#                     subject_to='fsaverage',  # Default: 'fsaverage'
#                     subjects_dir=subjects_dir,  # Default: None
#                     src=src_vol,  # Default: None
#                     niter_affine=(100, 100, 10),  # Default: (100, 100, 10)
#                     niter_sdr=(5, 5, 3),  # Default: (5, 5, 3)
#                     spacing=5)  # Default: 5

# Apply morph
# stc_fsaverage = morph(stc_vol)

# Make NIfTI volume variant 1
# img = morph(stc_vol,
#             as_volume=True,  # Default: False
#             mri_resolution=True,  # Default: False
#             mri_space=True, # Default: True
#             format='nifti2',  # Default: 'nifti1'
#             apply_morph=True)  # Default: False

# Make NIfTI volume variant 2
# img = morph.as_volume(stc_fsaverage,
#                       mri_resolution=(3., 3., 3.),  # iso voxel size 3 mm
#                       mri_space=True)

# Save morph to disk
# morph.save('my-favorite-morph.h5')

# Read morph from disk
# morph = read_source_morph('my-favorite-morph.h5')

# Shortcuts
# stc_fsaverage = SourceMorph(src=src_vol)(stc_vol)
# img = SourceMorph(src=src_vol)(stc_vol, as_volume=True, mri_resolution=True)

###############################################################################
# Setting up :class:`mne.SourceMorph`
# -----------------------------------
#
# SourceMorph is a class that computes a morph operation from one subject to
# another depending on the underlying data. The result will be an instance of
# :class:`mne.SourceMorph`, that contains the mapping between the two spaces.
# At the very least the source space corresponding to the source estimate that
# is going to be morphed, has to be provided. Since stored data from both
# subjects of reference will be used, it is necessary to ensure that
# subject_from and subject_to, as well as subjects_dir are correctly set.

# SourceMorph initialization If src is not provided, the morph will not be
# pre-computed but instead will be prepared for morphing when calling. This
# works only with (Vector)SourceEstimate

morph_surf = SourceMorph(subject_from='sample',
                         subject_to='fsaverage',
                         subjects_dir=subjects_dir)

# Ideally subject_from can be inferred from src, subject_to is 'fsaverage' by
# default and subjects_dir is set in the environment. In that case SourceMorph
# can be initialized taking only src as parameter.

morph_vol = SourceMorph(subject_from='sample',
                        subject_to='fsaverage',
                        subjects_dir=subjects_dir,
                        spacing=(3., 3., 3.),  # grid spacing (3., 3., 3.) mm
                        src=src_vol)

###############################################################################
# Applying an instance of :class:`mne.SourceMorph`
# ------------------------------------------------
#
# Once we computed the morph for our respective dataset, we can morph the data,
# by giving it as an argument to the SourceMorph instance. This operation
# applies pre-computed transforms to stc.

stc_surf_m = morph_surf(stc_surf)  # SourceEstimate | VectorSourceEstimate
stc_vol_m = morph_vol(stc_vol)  # VolSourceEstimate


###############################################################################
# Reading and writing :class:`mne.SourceMorph` from and to disk
# -------------------------------------------------------------
#
# An instance of SourceMorph can be saved, by calling
# :meth:`morph.save <mne.SourceMorph.save>`. This methods allows for
# specification of a filename. The morph will be save in ".h5" format. If no
# file extension is provided, "-morph.h5" will be appended to the respective
# defined filename.
# In turn, reading a saved source morph can be achieved by using
# :func:`mne.read_source_morph`.

# morph_vol.save('my-file-name')

# -morph.h5 was attached because no file extension was provided when saving
# morph_vol = read_source_morph('my-file-name-morph.h5')

###############################################################################
# Spacing parameter of :class:`mne.SourceMorph`
# ---------------------------------------------
#
# When morphing a surface source estimate, spacing can be an int or a list of
# two arrays. In the first case the data will be morphed to an icosahedral
# mesh, having a resolution of spacing (typically 5) using
# :func:`mne.grade_to_vertices`. In turn, when morphing a volumetric source
# estimate, spacing can be a tuple of float representing the voxel size in each
# spatial dimension or a single value (int or float) to represent the very same
# but assigning equal values to all spatial dimensions. Voxel size referring to
# the spacing of the reference volumes when computing the volumetric morph in
# mm. Note that voxel size is inverse related to computation time and accuracy
# of the morph operation.
# Changing the spacing for a volumetric morph estimation, does not affect the
# later resolution of the source estimate sfter applying the morph. It is
# rather the resolution of morph estimation and hence should increased when
# aiming for more precision. The default is an isometric voxel size of 5 mm.
# In general it might be advisable to use a spacing that is smaller or equal to
# the actual grid spacing of the source estimate.

# Estimate non-linear volumetric morph based on grid spacing of (7., 7., 7.) mm

# morph = SourceMorph(src=src_vol, spacing=(7., 7., 7.))  # equiv. to spacing=7

###############################################################################
# niter parameters of :class:`mne.SourceMorph`
# ---------------------------------------------
#
# Additionally, compuation time and accuracy of the respective volumetric morph
# result, depend on the number of iterations per step of optimization. Under
# the hood, an Affine transformation is computed based on the mutual
# information. This metric relates structural changes in image intensity
# values. Because different still brains expose high structural similarities
# this method works quite well to relate corresponding features [1]_. The
# nonlinear transformations will be performed as Symmetric Diffeomorphic
# Registration (sdr) using the cross-correlation metric [2]_.
# Both optimization procedures are performed in "levels", passing the result
# from the first level of refinement to the next and so on. For each level, the
# number of iterations to optimize the alignment, can be defined. This is done
# be assigning a tuple to niter_affine and niter_sdr. Each tuple contains as
# many values, as desired levels of refinement and each value, represents the
# number of iterations for the respective level. The default is
# niter_affine=(100, 100, 10) and niter_sdr=(5, 5, 3). Both algorithms will be
# performed using 3 levels of refinement each and the corresponding number of
# iterations.

# Estimate non-linear volumetric morph based on grid spacing of (7., 7., 7.) mm
# and a reduced number of iterations. Note the difference in computation time.

# morph = SourceMorph(src=src_vol,
#                     spacing=(7., 7., 7.),
#                     niter_affine=(10, 10, 10),  # 3 levels a 10 iterations
#                     niter_sdr=(3, 3))  # 2 levels a 3 iterations


###############################################################################
# Transforming :class:`mne.VolSourceEstimate` into NIfTI
# ------------------------------------------------------
#
# In case of the volume source estimate, we can further ask the morph to output
# a volume of our data in the new space. We do this by calling the
# :meth:`morph.as_volume <mne.SourceMorph.as_volume>`. Note, that un-morphed
# source estimates still can be converted into a NIfTI by using
# :meth:`stc.as_volume <mne.VolSourceEstimate.as_volume>`. The shape of the
# output volume can be modified by providing the argument mri_resolution. This
# argument can be boolean, a tuple or an int. If mri_resolution=True, the MRI
# resolution, that was stored in src will be used. Setting mri_resolution to
# False, will export the volume having voxel size corresponding to the spacing
# of the computed morph. Setting a tuple or single value, will cause the output
# volume to expose a voxel size of that values in mm.

# img_mri_res = morph_vol.as_volume(stc_vol_m, mri_resolution=True)

# img_morph_res = morph_vol.as_volume(stc_vol_m, mri_resolution=False)

# img_any_res = morph_vol.as_volume(stc_vol_m, mri_resolution=3)



###############################################################################
# Plot results
# ------------

# Plot morphed volume source estiamte

# Load fsaverage anatomical image
t1_fsaverage = nib.load(fname_t1_fsaverage)

# Initialize figure
fig, axes = plt.subplots()
fig.subplots_adjust(top=0.8, left=0.1, right=0.9, hspace=0.5)
fig.patch.set_facecolor('white')

# Setup nilearn plotting
display = plot_glass_brain(t1_fsaverage,
                           display_mode='ortho',
                           cut_coords=[0., 0., 0.],
                           draw_cross=False,
                           axes=axes,
                           figure=fig,
                           annotate=False)

# Transform into volume time series and use first one
overlay = index_img(morph_vol.as_volume(stc_vol_m, mri_resolution=True), 0)

display.add_overlay(overlay, alpha=0.75)
display.annotate(size=8)
axes.set_title('Morphed to fsaverage', color='white', fontsize=20)

plt.text(plt.xlim()[1], plt.ylim()[0], 't = 0.09s', color='white')
plt.show()

del stc_vol_m, morph_vol, morph_surf

# Plot morphed surface source estiamte

surfer_kwargs = dict(
    hemi='lh', subjects_dir=subjects_dir,
    clim=dict(kind='value', lims=[8, 12, 15]), views='lateral',
    initial_time=0.09, time_unit='s', size=(800, 800),
    smoothing_steps=5)
brain = stc_surf_m.plot(**surfer_kwargs)
brain.add_text(0.1, 0.9, 'Morphed to fsaverage', 'title', font_size=20)

del stc_surf_m

###############################################################################
# .. note:: If you are using an IIR filter, :func:`mne.filter.create_filter`
#           will not print a filter length and transition bandwidth to the log.
#           Instead, you can specify the roll off with the `iir_params`
#           argument or stay with the default, which is a 4th order
#           (Butterworth) filter.
#
# Passband ripple and stopband attenuation
# ----------------------------------------
#
# +-------------------------+-----------------+----------------------+
# | Name of window function | Passband ripple | Stopband attenuation |
# +=========================+=================+======================+
# | Hann                    | 0.0545 dB       | 44 dB                |
# +-------------------------+-----------------+----------------------+
# | Hamming                 | 0.0194 dB       | 53 dB                |
# +-------------------------+-----------------+----------------------+
# | Blackman                | 0.0017 dB       | 74 dB                |
# +-------------------------+-----------------+----------------------+
#
#
# Summary
# =======
#
#
# References
# ==========
#
# .. [1] Gramfort, A., Luessi, M., Larson, E., Engemann, D. A., Strohmeier, D.,
#        Brodbeck, C., ... & Hämäläinen, M. (2013). MEG and EEG data analysis
#        with MNE-Python. Frontiers in neuroscience, 7, 267.
# .. [2] Avants, B. B., Epstein, C. L., Grossman, M., & Gee, J. C. (2009).
#        Symmetric diffeomorphic image registration with cross-correlation:
#        evaluating automated labeling of elderly and neurodegenerative brain.
#        Medical image analysis, 12(1), 26-41.
# .. [3] Viola, P., & Wells III, W. M. (1997). Alignment by maximization of
#        mutual information. International journal of computer vision, 24(2),
#        137-154.
# .. [4] Mattes, D., Haynor, D. R., Vesselle, H., Lewellen, T. K., & Eubank, W.
#        (2003). PET-CT image registration in the chest using free-form
#        deformations. IEEE transactions on medical imaging, 22(1), 120-128.
# .. [5] Greve D. N., Van der Haegen L., Cai Q., Stufflebeam S., Sabuncu M.
#        R., Fischl B., Brysbaert M. A Surface-based Analysis of Language
#        Lateralization and Cortical Asymmetry. Journal of Cognitive
#        Neuroscience 25(9), 1477-1492, 2013.
# .. _dipy example: http://nipy.org/dipy/examples_built/syn_registration_3d.html  # noqa