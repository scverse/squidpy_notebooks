#!/usr/bin/env python
"""
Histogram features
------------------
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

import scanpy as sc
import squidpy as sq

###############################################################################
# Lets load a fluorescence Visium dataset and calculate bin-counts (3 bins) of channels 0 and 1.

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
# :ref:`sphx_glr_auto_tutorials_tutorial_napari.py` to learn how to use our interactive :mod:`napari` plugin.
# With these features we can e.g. appreciate the detailed distribution of
# intensity values of channel 0 (DAPI stain) on the different bins.

sc.pl.spatial(
    sq.pl.extract(adata, "histogram_features"),
    color=[None, "histogram_ch_0_bin_0", "histogram_ch_0_bin_1", "histogram_ch_0_bin_2"],
    bw=True,
)
