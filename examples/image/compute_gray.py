"""
Converting to Grayscale
------------------

This example shows how to use :func:`squidpy.im.process_img` to convert and image layer
of :class:`squidpy.im.ImageContainer` to grayscale.

We use the argument ``processing="gray"`` to smooth the image.
This calls :func:`skimage.color.rgb2gray` in the background.

See also :ref:`sphx_glr_auto_examples_image_compute_smooth.py`
and :ref:`sphx_glr_auto_examples_image_compute_process_hires.py`
"""

import os

import squidpy as sq

import matplotlib.pyplot as plt

# get spatial dataset including hires tissue image
# set path to dataset
BASE_PATH = "/Users/hannah.spitzer/projects/spatial_scanpy/data"
dataset_folder = os.path.join(BASE_PATH, "20191205_10XVisium_MouseBrainCoronal_giovanni.palla")
# load data
img = sq.im.ImageContainer(os.path.join(dataset_folder, "V1_Adult_Mouse_Brain_image.tif"))

# %%
# First, we crop a smaller image to smooth.
# This is only to speed things up, :func:`squidpy.im.process_img` can also process very large images
# (see :ref:`sphx_glr_auto_examples_image_compute_process_hires.py`.)
crop = img.crop_corner(4000, 4000, 250, 250)

# %%
# Convert to image to grayscale and plot the result.
# With the argument ``img_id`` we can select the image layer that should be processed.
# When converting to grayscale, the channel dimensions change from 3 to 1.
# By default, the name of the resulting channel dimension will be ``{{original_channel_name}}_gray``.
# Use the argument ``channel_id`` to set a new channel name explicitely.
# By default, the resulting image is saved in the layer ``image_gray`.
# This behaviour can be changed with the arguments ``copy`` and ``key_added``.

sq.im.process_img(crop, img_id="image", processing="gray")

fig, axes = plt.subplots(1, 2)
axes[0].imshow(crop["image"])
axes[1].imshow(crop["image_gray"], cmap="gray")
for ax in axes:
    ax.axis("off")
