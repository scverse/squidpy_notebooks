#!/usr/bin/env python
r"""
Extract texture features
------------------------

This example shows how to extract texture features from the tissue image.

Textures features give give a measure of how the image intensity at different distances and angles varies by
calculating a grey-level co-occurrence matrix (`GLCM <https://en.wikipedia.org/wiki/Co-occurrence_matrix>`_).
The GLCM includes the number of times that grey-level :math:`j` occurs at a distance :math:`d`
and at an angle :math:`\\theta` from grey-level :math:`i`.
From this data, different features (``props``) are calculated.
See also :func:`skimage.feature.greycomatrix`.

Texture features are calculated by using ``features = 'texture'``, which calls
:func:`squidpy.im.ImageContainer.features_texture`.
In addition to ``feature_name`` and ``channels``, we can also specify the following ``features_kwargs``:

    - ``distances`` - distances that are taken into account for finding repeating patterns.
    - ``angles`` - range on which values are binned. Default is the whole image range.
    - ``props`` - texture features that are extracted from the GLCM.

.. seealso::

    See :ref:`sphx_glr_auto_examples_image_compute_features.py` for general usage of
    :func:`squidpy.im.calculate_image_features`.
"""

import scanpy as sc
import squidpy as sq

###############################################################################
# Let's load the fluorescence Visium dataset and calculate texture features with default ``features_kwargs``.
#
# Note that for texture features it may make sense to compute them over a larger crop size to include more context,
# e.g., ``spot_scale = 2`` or ``spit_scale = 4`` which will extract crops with double or four times the radius
# than the original Visium spot size.
# For more details on the image cropping, see :ref:`sphx_glr_auto_examples_image_compute_crops.py`.

# get spatial dataset including high-resolution tissue image
img = sq.datasets.visium_fluo_image_crop()
adata = sq.datasets.visium_fluo_adata_crop()

# calculate texture features and save in key "texture_features"
sq.im.calculate_image_features(
    adata,
    img,
    features="texture",
    key_added="texture_features",
    spot_scale=2,
    show_progress_bar=False,
)
###############################################################################
# The result is stored in ``adata.obsm['texture_features']``.
adata.obsm["texture_features"].head()

###############################################################################
# Use :func:`squidpy.pl.extract` to plot the texture features on the tissue image or have a look at
# `our interactive visualization tutorial <../../external_tutorials/tutorial_napari.ipynb>`_ to learn
# how to use our interactive :mod:`napari` plugin.
# Here, we show the contrast feature for channels 0 and 1.
# The two stains, DAPI in channel 0, and GFAP in channel 1 show different regions of high contrast.
sc.pl.spatial(
    sq.pl.extract(adata, "texture_features"),
    color=[None, "texture_ch-0_contrast_dist-1_angle-0.00", "texture_ch-1_contrast_dist-1_angle-0.00"],
    bw=True,
)
