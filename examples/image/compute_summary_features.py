"""
Summary features
----------------

Here, we use :func:`squidpy.im.calculate_image_features` to extract summary features from the tissue image.
Have a look at :ref:`sphx_glr_auto_examples_image_compute_features.py`
for the general usage of :func:`squidpy.im.calculate_image_features`.

Summary features give a good overview over the intensity of each image channel at the location of the visium spots.
They are calculated by using ``features='summary'``.

In addition to ``feature_name`` and ``channels`` we can specify the following ``features_kwargs``:

- ``quantiles``: Quantiles that are computed. By default, the 0.9th, 0.5th, and 0.1th quantiles are calculated
- ``mean``: Compute mean. Off by default
- ``std``: Compute std deviation. Off by default.

"""

import scanpy as sc
import squidpy as sq

###############################################################################
# First, we load a fluorescence visisum dataset.

# get spatial dataset including high-resolution tissue image
img = sq.datasets.visium_fluo_image_crop()
adata = sq.datasets.visium_fluo_adata_crop()


###############################################################################
# Then, we calculate the 0.9th quantile and mean for the visium spots of the fluorescence channels 0 (DAPI)
# and 1 (GFAP).
# In order to only get statistics of the tissue underneath the spots, we use the argument ``mask_circle=True``.
# When not setting this flag, statistics are calculated using a square crop centered on the spot.

# calculate summary features and save in key "summary_features"
sq.im.calculate_image_features(
    adata,
    img,
    features="summary",
    features_kwargs={
        "summary": {
            "mean": True,
            "quantiles": [0.9],
            "channels": [0, 1],
        }
    },
    key_added="summary_features",
    mask_circle=True,
)

###############################################################################
# The result is stored in `adata.obsm['summary_features']`
adata.obsm["summary_features"].head()

###############################################################################
# Use :func:`squidpy.pl.extract` to plot the summary features on the tissue image or have a look at
# :ref:`sphx_glr_auto_tutorials_tutorial_napari.py` to learn how to use our interactive napari plugin.
# Note, how the spatial distribution of channel means is different for fluorescence channels 0 (DAPI stain)
# and 1 (GFAP stain).

sc.pl.spatial(sq.pl.extract(adata, "summary_features"), color=[None, "summary_mean_ch_0", "summary_mean_ch_1"], bw=True)
