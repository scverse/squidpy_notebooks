"""
Smoothing an Image
------------------

This example shows how to use :func:`squidpy.im.process_img` to smooth an image layer of
:class:`squidpy.im.ImageContainer`.

We use the argument ``processing="smooth"`` to smooth the image.
This calls :func:`skimage.filters.gaussian` in the background.
By providing the argument ``processing_kwargs`` we can pass keyword arguments to the wrapped function.
This allows us to set the width of the gaussian kernel, ``sigma``, used for smoothing.

See also :ref:`sphx_glr_auto_examples_image_compute_gray.py` and
:ref:`sphx_glr_auto_examples_image_compute_process_hires.py`
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
# Smooth the image with ``"sigma" = 2``.
# With the argument ``img_id`` we can select the image layer that should be processed.
# By default, the resulting image is saved in the layer ``image_smooth`.
# This behaviour can be changed with the arguments ``copy`` and ``key_added``.

sq.im.process_img(crop, img_id="image", processing="smooth", processing_kwargs={"sigma": 2})

# %%
# Now we can plot the result
fig, axes = plt.subplots(1, 2)
axes[0].imshow(crop["image"])
axes[1].imshow(crop["image_smooth"])
for ax in axes:
    ax.axis("off")
