# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: sphinx
#       format_version: '1.1'
#       jupytext_version: 1.9.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

"""
Extract image features using a pretrained deep learning model
=============================================================

In recent years, deep learning models have become increasing popular for computer vision tasks.
It is possible to use a pretraine neural network to extract features from an image.
These DL-powered features can be very powerful and provide a semantically meaningful compression of the image.

Here, we show how to use the feature extraction pipeline defined in :mod:`squidpy.im` to extract features from
a pretrained neural network.
For this, we use a ResNet50 ([ResNet16]_) with weights pretrained on `ImageNet`_.

To execute this notebook, you need to install ``tensorflow`` in addition to ``squidpy``.
TODO more info on how to install

.. _ImageNet: http://www.image-net.org
"""

import time

import tensorflow as tf

import scanpy as sc
import anndata as ad
import squidpy as sq

import pandas as pd

###############################################################################
# Load data and prepare ResNet
# ----------------------------
#
# Load a mouse brain visium dataset with an H&E stained tissue image.

img = sq.datasets.visium_hne_image_crop()
adata = sq.datasets.visium_hne_adata_crop()

###############################################################################
# By default, :func:`squidpy.im.calculate_image_features` uses image crops centered on the spots that
# have the same diameter as the spots.
# This behaviour is configureable by specifying the arguments ``size`` and ``scale``.
# Here, we are going to use slightly lower resolution (``scale=0.5``) image crops
# that are double the size of the spots (``size=2``) for more context.
# For more information see also :ref:`sphx_glr_auto_examples_image_compute_features.py`.
#
# To determine the shape of the image crops that will be passed to the DL feature extractor, we use
# :func:`squidpy.im.generate_spot_crops` with the desired parameters.

size = 2
scale = 0.5
crop_shape = next(img.generate_spot_crops(adata, size=size, scale=scale))[0]["image"].shape

print(f"the image crop shape is {crop_shape}")

###############################################################################
# Lets load the `ResNet50` with pretrained weights.
# We pass the previously determined ``input_shape`` and specify ``include_top=False``,
# as we are going to use the outputs of the last convolutional layer as image features.

resnet = tf.keras.applications.ResNet50(include_top=False, weights="imagenet", input_shape=crop_shape)


###############################################################################
# ``resnet`` is a keras model. To make predictions, simply pass an image to it.
# The output of this will be a 4 dimensional tensor that we will flatten to get a 1d feature vector.
#
# Let us define ``resnet_features()`` as a helper function to extract the features from a given array using resnet.
# In the following, we will use this function with :func:`squidpy.im.calculate_image_features` to extract features
# from each spot location.


def resnet_features(arr):
    """Extract features using resnet."""
    import numpy as np

    # add batch dimension to arr
    arr = arr[np.newaxis, :, :, :]
    # predict features
    features = resnet(arr)
    # to numpy ans flatten
    features = features.numpy().flatten()
    return features


###############################################################################
# Extract ResNet features
# -----------------------
#
# Now we are ready to extract features.
# We are using the `"custom"` feature extractor, that allows calculation of features with a custom
# feature extraction function.
#
# On a CPU, executing the following will take roughly 1.5mins.

start_time = time.time()
sq.im.calculate_image_features(
    adata, img, features="custom", features_kwargs={"custom": {"feature_fn": resnet_features}}, size=size, scale=scale
)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.2f}s")


###############################################################################
# Visualise and compare feature and gene clusters
# -----------------------------------------------
#
# Define ``cluster_features()`` as a helper function for clustering a feature matrix.


def cluster_features(features: pd.DataFrame, like=None):
    """Calculate leiden clustering of features.

    Specify filter of features using `like`.
    """
    # filter features
    if like is not None:
        features = features.filter(like=like)
    # create temporary adata to calculate the clustering
    adata = ad.AnnData(features)
    # adata.var_names_make_unique()
    # important - feature values are not scaled, so need to scale them before PCA
    sc.pp.scale(adata)
    # calculate leiden clustering
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)

    return adata.obs["leiden"]


###############################################################################
# Cluster ResNet features and compare with a gene-based clustering

adata.obs["resnet_cluster"] = cluster_features(adata.obsm["img_features"])
sc.set_figure_params(facecolor="white")
sc.pl.spatial(adata, color=["resnet_cluster", "cluster"])

###############################################################################
# The clustering based on ResNet features respects uniform areas on the H&E stain.
# For example the pyramidial layer of the hippocampus is assigned to one cluster.
