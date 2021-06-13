#!/usr/bin/env python
"""
Extract custom features
-----------------------

This example shows how to extract features from the tissue image using a custom function.

The custom feature calculation function can be any python function that takes an image as input, and
returns a list of features.
Here, we show a simple example by defining a function to calculate the mean of the images.

Custom features are calculated by using ``features = 'custom'``, which calls
:func:`squidpy.im.ImageContainer.features_custom`.
In addition to ``feature_name`` and ``channels`` we can specify the following ``features_kwargs``:

    - ``func`` - custom feature extraction function.
    - ``additional_layers`` - names of image layers that should be passed to ``func`` together with ``layer``.
    - additional keyword arguments for ``func``.

.. seealso::

    See :ref:`sphx_glr_auto_examples_image_compute_features.py` for general usage of
    :func:`squidpy.im.calculate_image_features`.
"""

import scanpy as sc
import squidpy as sq

###############################################################################
# Let's load the H&E Visium dataset.

# get spatial dataset including high-resolution tissue image
img = sq.datasets.visium_hne_image_crop()
adata = sq.datasets.visium_hne_adata_crop()


###############################################################################
# Define a custom feature extraction function.
def mean_fn(arr):
    """Compute mean of arr."""
    import numpy as np

    return np.mean(arr)


###############################################################################
# Now we can extract features using `mean_fn` by providing it within ``features_kwargs``.
sq.im.calculate_image_features(
    adata,
    img,
    features="custom",
    features_kwargs={"custom": {"func": mean_fn}},
    key_added="custom_features",
    show_progress_bar=False,
)

###############################################################################
# The result is stored in ``adata.obsm['custom_features']``.
adata.obsm["custom_features"].head()

###############################################################################
# Use :func:`squidpy.pl.extract` to plot the histogram features on the tissue image or have a look at
# `our interactive visualization tutorial <../../external_tutorials/tutorial_napari.ipynb>`_ to learn
# how to use our interactive :mod:`napari` plugin.
sc.pl.spatial(
    sq.pl.extract(adata, "custom_features"),
    color=[None, "mean_fn_0"],
    bw=True,
)


###############################################################################
# You can also pass more than one image layer to the custom feature extraction function.
# For this, specify the necessary additional layer names using ``additional_layers`` in ``features_kwargs``.
# The specified image layers will be passed to the custom feature extraction function.
#
# Here, we show this behavior by defining a feature extraction function that sums two image layers:
def sum_fn(arr, extra_layer):
    """Compute sum of two image layers."""
    import numpy as np

    return np.sum(arr + extra_layer)


img.add_img(img["image"].values, layer="extra_layer")

sq.im.calculate_image_features(
    adata,
    img,
    layer="image",
    features="custom",
    features_kwargs={"custom": {"func": sum_fn, "additional_layers": ["extra_layer"]}},
    key_added="custom_features",
    show_progress_bar=False,
)
