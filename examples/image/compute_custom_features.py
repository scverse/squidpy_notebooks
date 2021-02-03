"""
Custom features
---------------

Here, we use :func:`squidpy.im.calculate_image_features` to extract features from the tissue image
using a custom function.
Have a look at :ref:`sphx_glr_auto_examples_image_compute_features.py` for the general usage of
:func:`squidpy.im.calculate_image_features`.

The custom feature calculation function can be any python function that takes an image as input, and
returns a list of features.

Here, we show a simple example by defining a function to calculate the mean of the images.
"""

import scanpy as sc
import squidpy as sq

###############################################################################
# Lets load a H&E visisum dataset

# get spatial dataset including high-resolution tissue image
img = sq.datasets.visium_hne_image_crop()
adata = sq.datasets.visium_hne_adata_crop()

###############################################################################
# Define a custom feature extraction function


def feature_fn(arr):
    """Compute mean of arr."""
    import numpy as np

    return np.mean(arr)


###############################################################################
# Now we can extract features using ``feature_fn``

sq.im.calculate_image_features(
    adata, img, features="custom", features_kwargs={"custom": {"feature_fn": feature_fn}}, key_added="custom_features"
)

###############################################################################
# The result is stored in ``adata.obsm['custom_features']``.

adata.obsm["custom_features"].head()

###############################################################################
# Use :func:`squidpy.pl.extract` to plot the histogram features on the tissue image or have a look at
# :ref:`sphx_glr_auto_tutorials_tutorial_napari.py` to learn how to use our interactive napari plugin.

sc.pl.spatial(
    sq.pl.extract(adata, "custom_features"),
    color=[None, "custom_0"],
    bw=True,
)
