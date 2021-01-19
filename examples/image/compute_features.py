"""
Extract Image Features
----------------------
In addition to the spatial gene expression values, visium datasets also contain high-resolution images of the tissue.
In this notebook we extract features from this image using :func:`squidpy.im.calculate_image_features`
and create a ``obs x features`` matrix that can be analysed together with
the ``obs x genes`` spatial gene expression matrix.

To compute features for each visium spot (``obs``), we extract image crops from the tissue image centered on each spot.
When extracting the crops, we can specify the size and scale of the crops and optionally mask a circle
to ensure that only tissue underneath the round visium spots is taken into account to compute the features.
See also :ref:`sphx_glr_auto_examples_image_compute_crops.py`.

The extracted crops are then used to compute features.
We provide different feature extractors that are described in more detail in the following examples:
- summary statistics of each color channel
  (:ref:`sphx_glr_auto_examples_image_compute_summary_features.py`)
- texture features based on repeating patterns
  (:ref:`sphx_glr_auto_examples_image_compute_texture_features.py`)
- color histogram features using counts in bins of each channel's histogram
  (:ref:`sphx_glr_auto_examples_image_compute_histogram_features.py`)
- number and size of objects from a binary segmentation layer
  (:ref:`sphx_glr_auto_examples_image_compute_segmentation_features.py`)

"""
# import modules
import os

import squidpy as sq

import scanpy as sc

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

# get spatial dataset including hires tissue image
img = sq.im.ImageContainer(os.path.expanduser("~/.cache/squidpy/tutorial_data/visium_hne_crop.tiff"))
adata = sc.read(os.path.expanduser("~/.cache/squidpy/tutorial_data/visium_hne_crop.h5ad"))

# %%
# The high resolution tissue image is contained in ``img['image']``,
# and the spot locations in tissue image pixel-space are located in ``adata.obsm['spatial']``.
# We can plot the spots overlayed on a lower-resolution version of the tissue image contained in adata.

np.set_printoptions(threshold=10)
print(img)
print(adata.obsm["spatial"])

sc.pl.spatial(adata, add_outline=True)

# %%
# Using this information, we can now extract features from the tissue underneath each spot by calling
# :func:`squidpy.im.calculate_image_features`.
# This function takes both `adata` and `img` as input, and will write the resulting ``obs x features`` matrix to
# ``adata.obsm[key]``.
# It contains several arguments to modify its behaviour.
# With these arguments you can
# - specify the image used for feature calculation using ``img_id``,
# - specify the type of features that should be calculated using ``features`` and ``features_kwargs``,
# - specify how the crops used for feature calculation look like using ``kwargs``,
# - specify parallelization options using ``n_jobs``, ``backend``, ``show_progress_bar``, and
# - specify how the data that is returned using ``key_added``, ``copy``.
#
# Let us first calculate summary features and save the result in ``adata.obsm['features']``.

sq.im.calculate_image_features(adata, img, features="summary", key_added="features")

# show the calculated features
print(f"calculated features: {list(adata.obsm['features'].columns)}")
print(adata.obsm["features"].head())

# %%
# To visualize the features, we can use :func:`squidpy.pl.extract` to plot the texture features on the tissue image.
# See :ref:`sphx_glr_auto_examples_plotting_compute_extract.py` for more details on this function.
#
# Here, we plot the median value of channel 0 (``summary_quantile_0.5_ch_0``).

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.pl.spatial(sq.pl.extract(adata, "features"), color=[None, "summary_quantile_0.5_ch_0"])

# %%
# ## Speed up feature extraction
# Speeding up the feature extraction is easy.
# Just set the ``n_jobs`` flag to the number of jobs that should be used by :func:`squidpy.im.calculate_image_features`.

# extract features by using 4 jobs
sq.im.calculate_image_features(adata, img, features="summary", key_added="features", n_jobs=4)

# %%
# ## Specify crop appearance
# Features are extracted from image crops that are centered on the visium spots
# (see also :ref:`sphx_glr_auto_examples_image_compute_crops.py`).
# By default, the crops have the same size as the spot, are not scaled and not masked.
# We can use the `mask_circle`, `scale`, and `size` arguments to change how the crops are generated.
#
# - Use `mask_circle=True, scale=1, size=1`, if you would like to get features that are calculated only from tissue
#   in a visium spot
# - Use `scale=X`, with `X < 1`, if you would like to downscale the crop before extracting the features
# - Use `size=X`, with `X > 1`, if you would like to extract crops that are X-times the size of the visium spot
#
# Let us extract masked and scaled features and compare them

# We subset adata to the first 50 spots to make the computation of features fast.
# Skip this step if you want to calculate features from all spots
adata_sml = adata[0:50].copy()

# calculate default features
sq.im.calculate_image_features(adata_sml, img, features=["summary", "texture", "histogram"], key_added="features")
# calculate features with masking
sq.im.calculate_image_features(
    adata_sml, img, features=["summary", "texture", "histogram"], key_added="features_masked", mask_circle=True
)
# calculate features with scaling and larger context
sq.im.calculate_image_features(
    adata_sml,
    img,
    features=["summary", "texture", "histogram"],
    key_added="features_scaled",
    mask_circle=True,
    size=2,
    scale=0.5,
)

# plot distribution of median for different cropping options
sns.distplot(adata_sml.obsm["features"]["summary_quantile_0.5_ch_0"], label="features")
sns.distplot(adata_sml.obsm["features_masked"]["summary_quantile_0.5_ch_0"], label="features_masked")
sns.distplot(adata_sml.obsm["features_scaled"]["summary_quantile_0.5_ch_0"], label="features_scaled")
_ = plt.legend()

# %%
# The masked features have lower median values, because the area outside the circle is masked with zeros.
