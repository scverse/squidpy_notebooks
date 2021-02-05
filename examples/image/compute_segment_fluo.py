"""
Cell-segmentation
------------------

We can use the high resolution tissue images to segment nuclei.
This information can be used to compute additional image features like cell count and cell size per spot
(see :ref:`sphx_glr_auto_examples_image_compute_segmentation_features.py`).
This example shows how to use :func:`squidpy.im.segment_img` and explains the parameters you can use.

We provide two segmentation models :class:`squidpy.im.SegmentationModelBlob`
and :class:`squidpy.im.SegmentationModelWatershed`.
In addition, you can use your own pre-trained :mod:`tensorflow.keras` model to perform the segmentation
utilising :class:`squidpy.im.SegmentationModelTensorflow`.

Note that when using the provided segmentation models ``"blob"`` and ``"watershed"``, the quality of the
cell-segmentation depends on the quality of your tissue images.
In this example we use the DAPI stain of a fluorescence dataset that clearly shows the nuclei to do the segmentation.
For harder cases, you may want to provide your own pre-trained segmentation model.

See :ref:`sphx_glr_auto_examples_image_compute_segment_hne.py` for an example of how to
calculate a cell-segmentation of an H&E stain.
"""

import squidpy as sq

import numpy as np

import matplotlib.pyplot as plt

# load fluorescence tissue image
img = sq.datasets.visium_fluo_image_crop()


###############################################################################
# We crop the image to a smaller segment.
# This is only to speed things up, :func:`squidpy.im.segment_img` can also process very large images
# (see :ref:`sphx_glr_auto_examples_image_compute_process_hires.py`.)
crop = img.crop_corner(1000, 1000, size=1000)

###############################################################################
# The tissue image in this dataset contains four fluorescence stains.
# The first one is DAPI, which we will use for the nuclei-segmentation.

fig, axes = plt.subplots(1, 3, figsize=(10, 20))
for i, ax in enumerate(axes):
    ax.imshow(crop["image"][:, :, i])
    ax.axis("off")

###############################################################################
# We segment the image with :func:`squidpy.im.segment` using watershed segmentation
# (``method="watershed"``).
# With the arguments ``image_id`` and ``channel`` we define the image layer and
# channel of this image that should be segmented.
#
# With ``kwargs`` we can provide keyword arguments to the segmentation model.
# For watershed segmentation, we need to set a threshold to create the mask image.
# You can either set a manual threshold, or use automated
# `Otsu thresholding <https://en.wikipedia.org/wiki/Otsu%27s_method>`_.
# For this florescence image example, Otsu thresholding works very well,
# thus we will use ``thresh=None``.
# See :ref:`sphx_glr_auto_examples_image_compute_segment_hne.py`
# for an example where we use a manually defined threshold.
#
# In addition, we can specify if the values greater or equal than
# the threshold should be in the mask (default)
# or if the values smaller to the threshold should be in the mask (``geq=False``).
#
sq.im.segment(img=crop, img_id="image", channel=0, method="watershed", thresh=None, geq=True)

###############################################################################
# The segmented crop is saved in the layer ``segmented_watershed``.
# This behaviour can be changed with the arguments ``copy`` and ``key_added``.
# The result of the segmentation is a label image that can be used to extract features like the
# number of cells from the image.
print(crop)
print(f"number of segments in crop: {len(np.unique(crop['segmented_watershed']))}")

fig, axes = plt.subplots(1, 2)
axes[0].imshow(crop["image"][:, :, 0])
axes[0].set_title("DAPI")
axes[1].imshow(crop["segmented_watershed"].squeeze(), cmap="jet", interpolation="none")
axes[1].set_title("segmentation")
for ax in axes:
    ax.axis("off")
