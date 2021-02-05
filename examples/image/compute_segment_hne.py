"""
Advanced Cell-segmentation for H&E stains
-----------------------------------------

This example shows how to use processing and segmentation functions to segment images with H&E stains.
For a general example of how to use :func:`squidpy.im.segment_img`
see :ref:`sphx_glr_auto_examples_image_compute_segment_fluo.py`.

Here, we attempt to segment a noisy H&E stain.
Note that we only provide very basic segmentation models.
If you require precise cell-segmentation and cell-counts, you might want to add more pre-processing
and / or use a pre-trained model to do the segmentation (using :class:`squidpy.im.SegmentationModelTensorflow`).
"""

import squidpy as sq

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

# load H&E stained tissue image and crop to a smaller segment
img = sq.datasets.visium_hne_image_crop()
crop = img.crop_corner(0, 0, size=1000)

###############################################################################
# Before segmenting the image, we do some preprocessing using :func:`squidpy.im.process`.

# smooth image
sq.im.process(crop, img_id="image", method="smooth", sigma=4)

# plot the result
fig, axes = plt.subplots(1, 2)
for img_id, ax in zip(["image", "image_smooth"], axes):
    ax.imshow(np.squeeze(crop[img_id]))
    ax.set_title(img_id)
    ax.axis("off")

###############################################################################
# We will use channel 0 to do the segmentation, as this channel contains most of
# the nuclei information within an H&E stain.
# Instead of using Otsh thresholding, we will define a manual fixed threshold.
# Note that using Otsu's method to determine the threshold also yields good results.
#
# Judging by the plot showing values smaller than 0.28, this threshold seems to be a good
# choice for this example.
fig, axes = plt.subplots(1, 3, figsize=(12, 5))
axes[0].imshow(crop["image_smooth"][:, :, 0], cmap="gray")
axes[0].axis("off")
axes[1].imshow(crop["image_smooth"][:, :, 0] < 0.36)
axes[1].axis("off")
_ = sns.histplot(np.array(crop["image_smooth"]).flatten(), bins=50, ax=axes[2])


###############################################################################
# We use :func:`squidpy.im.segment` with ``method="watershed"`` to do the segmentation.
# Since, opposite to the fluorescence DAPI stain, in the H&E stain, nuclei appear darker,
# we need to indicate the model that it should treat lower-intensity values as foreground.
# We do this by specifying the ``geq=False`` in the ``kwargs``.
sq.im.segment(img=crop, img_id="image_gray_smooth", method="watershed", thresh=0.28, geq=False)

###############################################################################
# The segmented crop is saved in the layer `segmented_watershed`.
# This behaviour can be changed with the arguments ``copy`` and ``key_added``.
# The result of the segmentation is a label image that can be used to extract features
# like the number of cells from the image.
print(crop)
print(f"number of segments in crop: {len(np.unique(crop['segmented_watershed']))}")

fig, axes = plt.subplots(1, 2)
axes[0].imshow(crop["image_gray_smooth"][:, :, 0])
axes[0].set_title("H&E")
axes[1].imshow(crop["segmented_watershed"].squeeze(), cmap="jet", interpolation="none")
axes[1].set_title("segmentation")
for ax in axes:
    ax.axis("off")
