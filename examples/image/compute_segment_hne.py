# %%
"""
Advanced Cell-segmentation for H&E stains
------------------

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
crop = img.crop_corner(0, 0, 1000, 1000)

# %%
# Before segmenting the image, we do some preprocessing using :func:`squidpy.im.process_img`.

# convert to grayscale
sq.im.process_img(crop, img_id="image", processing="gray")
# smooth image
sq.im.process_img(crop, img_id="image_gray", processing="smooth", sigma=4)

# plot the result
fig, axes = plt.subplots(1, 3, figsize=(12, 4)
for img_id, ax in zip(["image", "image_gray", "image_gray_smooth"], axes):
    ax.imshow(np.squeeze(crop[img_id]))
    ax.set_title(img_id)
    ax.axis("off")

# %%
# Finding a good threshold for the segmentation is more difficult than for a DAPI stain,
# as there is no distinct peak in the histogram. 
# Judging by the plot showing values smaller than 0.28, this threshold seems to be a good 
# choice for this example.
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
axes[0].imshow(crop["image_gray_smooth"][:, :, 0] < 0.28)
axes[0].axis("off")
_ = sns.histplot(np.array(crop["image_gray_smooth"]).flatten(), bins=50, ax=axes[1])


# %%
# We use :func:`squidpy.im.segment_img` with ``mode="watershed"`` to do the segmentation.
# Since, opposite to the fluorescence DAPI stain, in the H&E stain, nuclei appear darker,
# we need to indicate the model that it should treat lower-intensity values as foreground.
# We do this by specifying the ``geq=False`` in the ``kwargs``.
sq.im.segment_img(img=crop, img_id="image_gray_smooth", model_group="watershed", thresh=0.28, geq=False)

# %%
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
