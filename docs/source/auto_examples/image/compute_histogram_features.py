#!/usr/bin/env python
"""
Extract histogram features
--------------------------

This example shows how to extract histogram features from tissue image.

Histogram features give a more detailed view than summary features
(:ref:`sphx_glr_auto_examples_image_compute_summary_features.py`)
by computing a histogram of each image channel and returning bin-counts for each Visium spot.

In addition to ``feature_name`` and ``channels`` we can specify the following ``features_kwargs``:

    - ``bins`` - number of bins of the histogram, default is 10.
    - ``v_range`` - range on which values are binned, default is the whole image range.

.. seealso::

    See :ref:`sphx_glr_auto_examples_image_compute_features.py` for general usage of
    :func:`squidpy.im.calculate_image_features`.
"""

import squidpy as sq

###############################################################################
# Lets load the fluorescence Visium dataset and calculate bin-counts (3 bins) of channels 0 and 1.

# get spatial dataset including high-resolution tissue image
img = sq.datasets.visium_fluo_image_crop()
adata = sq.datasets.visium_fluo_adata_crop()

# calculate histogram features and save in key "histogram_features"
sq.im.calculate_image_features(
    adata,
    img,
    features="histogram",
    features_kwargs={"histogram": {"bins": 3, "channels": [0, 1]}},
    key_added="histogram_features",
)

###############################################################################
# The result is stored in ``adata.obsm['histogram_features']``.
adata.obsm["histogram_features"].head()

###############################################################################
# Use :func:`squidpy.pl.extract` to plot the histogram features on the tissue image or have a look at
# `our interactive visualisation tutorial <../../external_tutorials/tutorial_napari.ipynb>`_ to
# learn how to use our interactive :mod:`napari` plugin.
# With these features we can e.g. appreciate the detailed distribution of
# intensity values of channel 0 (DAPI stain) on the different bins.
sq.pl.spatial_scatter(
    sq.pl.extract(adata, "histogram_features"),
    color=[None, "histogram_ch-0_bin-0", "histogram_ch-0_bin-1", "histogram_ch-0_bin-2"],
    img_cmap="gray",
)
