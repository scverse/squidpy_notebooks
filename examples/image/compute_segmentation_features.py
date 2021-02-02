"""
Segmentation features
---------------------

Here, we use :func:`squidpy.im.calculate_image_features` to extract segmentation features from the tissue image.
Please have a look at :ref:`sphx_glr_auto_examples_image_compute_features.py` for the general usage of
:func:`squidpy.im.calculate_image_features`.

Segmentation features are computed using a label image where each object (nucleous)
in the image has a different number.
See :ref:`sphx_glr_auto_examples_image_compute_segment_fluo.py` for how to calculate such a label image

Features extracted from a nucleous segmentation range from the number of nuclei per image,
over numclei shapes and sizes, to the intensity of the input channels within the segmented objects.
They are very interpretable features and provide valuable additional information.
Use ``features='segmentation'`` to calculate the features.

In addition to ``feature_name`` and ``channels`` we can specify the following ``features_kwargs``:

- ``label_img_id``: name of label image layer in ``img``
- ``props``: Segmentation features that are calculated. See `properties` in :func:`skimage.measure.regionprops_table`.
- ``mean``: Return the mean of feature values per obs
- ``std``: Return the std of feature values per obs

"""
# sphinx_gallery_thumbnail_number = 2

import scanpy as sc
import squidpy as sq

import matplotlib.pyplot as plt

###############################################################################
# Lets load a fluorescence visium dataset

img = sq.datasets.visium_fluo_image_crop()
adata = sq.datasets.visium_fluo_adata_crop()

###############################################################################
# Before calculating segmentation features, we need to first calculate a segmentation
# using :func:`squidpy.im.segment_img`.

sq.im.segment_img(
    img=img, img_id="image", key_added="segmented_watershed", model_group="watershed", channel_idx=0, thresh=50000
)

###############################################################################
# Now we can calculate segmentation features.
# Here, we will calculate the following features:
# - number of nuclei (``label``)
# - mean area of nuclei (``area``)
# - mean intensity of channels 1 (anti-NEUN) and 2 (anti-GFAP) within nuclei (``mean_intensity``)
#
# We use ``mask_cicle=True`` to ensure that we are only extracting features from the tissue underneath each visium spot.
# For more details on the image cropping, see :ref:`sphx_glr_auto_examples_image_compute_crops.py`.

sq.im.calculate_image_features(
    adata,
    img,
    features="segmentation",
    key_added="segmentation_features",
    features_kwargs={
        "segmentation": {
            "label_img_id": "segmented_watershed",
            "props": ["label", "area", "mean_intensity"],
            "channels": [1, 2],
        }
    },
    mask_circle=True,
)

###############################################################################
# The result is stored in ``adata.obsm['segmentation_features']``

adata.obsm["segmentation_features"].head()

###############################################################################
# Use :func:`squidpy.pl.extract` to plot the texture features on the tissue image or have a look at
# :ref:`sphx_glr_auto_tutorials_tutorial_napari.py` to learn how to use our interactive napari plugin.
# Here, we show all calculated segmentation features.

# show all channels
fig, axes = plt.subplots(1, 3, figsize=(8, 4))
for i, ax in enumerate(axes):
    ax.imshow(adata.uns["spatial"]["V1_Adult_Mouse_Brain_Coronal_Section_2"]["images"]["hires"][:, :, i])
    ax.set_title(f"ch{i}")
    ax.axis("off")

# plot segmentation features
sc.pl.spatial(
    sq.pl.extract(adata, "segmentation_features"),
    color=[
        "segmentation_label",
        "segmentation_area_mean",
        "segmentation_mean_intensity_ch1_mean",
        "segmentation_mean_intensity_ch2_mean",
    ],
    bw=True,
    ncols=2,
)

###############################################################################
# `segmentation_label` shows the number of nuclei per spot and `segmentation_area_mean` the mean are of nuclei per spot.
# The remaining two plots show the mean intensity of channels 1 and 2 per spot.
# As the stains for channels 1 and 2 are specific to Neurons and Glial cells, respectively,
# these features show us Neuron and Glial cell dense areas.
