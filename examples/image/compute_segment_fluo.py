# %%
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

import seaborn as sns
import matplotlib.pyplot as plt

# load fluorescence tissue image
img = sq.datasets.visium_fluo_image_crop()


# %%
# We crop the image to a smaller segment.
# This is only to speed things up, :func:`squidpy.im.segment_img` can also process very large images
# (see :ref:`sphx_glr_auto_examples_image_compute_process_hires.py`.)
crop = img.crop_corner(1000, 1000, 1000, 1000)

# %%
# The tissue image in this dataset contains four fluorescence stains.
# The first one is DAPI, which we will use for the nuclei-segmentation.

fig, axes = plt.subplots(1, 3, figsize=(10, 20))
for i, ax in enumerate(axes):
    ax.imshow(crop["image"][:, :, i])
    ax.axis("off")

# %%
# For watershed segmentation, we need to set a threshold to create the mask image.
# The threshold should be chosen in such a way, that all nuclei are contained in the mask image.
_ = sns.histplot(np.asarray(crop["image"][:, :, 0]).flatten(), bins=50)

# %%
# The histogram of DAPI values shows a small peak at 60000 containing the nuclei.
# So, let us choose 50000 as a threshold for the segmentation method.
#
# We segment the image using the chosen threshold with :func:`squidpy.im.segment_img`.
# The argument ``image_id`` sets the image layer of img that should be segmented.
# Since we are segmenting the first channel, we will set ``channel_idx=0``.
# With the argument ``model_group`` we specify the model that we'd like to use for the segmentation.
# In our case this is ``"watershed"``.
# With ``kwargs`` we can provide keyword arguments to the segmentation model.
# For watershed, we need to set the threshold, ``thresh=50000``, as determined above.
# In addition, we can specify if the values greater or equal than the threshold should be in the mask (default)
# or if the values smaller to the threshold should be in the mask (``geq=False``).
sq.im.segment_img(img=crop, img_id="image", model_group="watershed", channel_idx=0, thresh=50000)

# %%
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
