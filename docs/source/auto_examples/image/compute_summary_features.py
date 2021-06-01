#!/usr/bin/env python
"""
Extract summary features
------------------------

This example shows how to extract summary features from the tissue image.

Summary features give a good overview over the intensity of each image channels at the location of the Visium spots.
They are calculated by using ``features = 'summary'``, which calls :func:`squidpy.im.ImageContainer.features_summary`.

In addition to ``feature_name`` and ``channels`` we can specify the following ``features_kwargs``:

- ``quantiles`` - quantiles that are computed. By default, the 0.9th, 0.5th, and 0.1th quantiles are calculated.

.. seealso::

    See :ref:`sphx_glr_auto_examples_image_compute_features.py` for general usage of
    :func:`squidpy.im.calculate_image_features`.
"""

import scanpy as sc
import squidpy as sq

###############################################################################
# First, we load a fluorescence Visium dataset.

# get spatial dataset including hires tissue image
img = sq.datasets.visium_fluo_image_crop()
adata = sq.datasets.visium_fluo_adata_crop()


###############################################################################
# Then, we calculate the 0.1th quantile, mean and standard deviation for the Visium spots
# of the fluorescence channels 0 (DAPI) and 1 (GFAP).
# In order to get statistics of only the tissue underneath the spots, we use the argument ``mask_circle = True``.
# When not setting this flag, statistics are calculated using a square crop centered on the spot.

# calculate summary features and save in key "summary_features"
sq.im.calculate_image_features(
    adata,
    img,
    features="summary",
    features_kwargs={
        "summary": {
            "quantiles": [
                0.1,
            ],
            "channels": [0, 1],
        }
    },
    key_added="summary_features",
    mask_circle=True,
    show_progress_bar=False,
)

###############################################################################
# The result is stored in ``adata.obsm['summary_features']``

adata.obsm["summary_features"].head()

###############################################################################
# Use :func:`squidpy.pl.extract` to plot the summary features on the tissue image or have a look at
# `our interactive visualisation tutorial <../../external_tutorials/tutorial_napari.ipynb>`_ to learn
# how to use our interactive :mod:`napari` plugin.
# Note how the spatial distribution of channel means is different for fluorescence channels 0 (DAPI stain)
# and 1 (GFAP stain).

sc.pl.spatial(sq.pl.extract(adata, "summary_features"), color=[None, "summary_ch-0_mean", "summary_ch-1_mean"], bw=True)
