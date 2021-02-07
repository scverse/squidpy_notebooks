#!/usr/bin/env python
"""
Convert to grayscale
--------------------
This example shows how to use :func:`sqduipy.im.process_img` to convert an image layer of
:class:`squidpy.im.ImageContainer` to grayscale.

We use the argument ``processing = 'gray'`` to convert the image. This calls :func:`skimage.color.rgb2gray`
in the background.

.. seealso::

    - :ref:`sphx_glr_auto_examples_image_compute_smooth.py`
    - :ref:`sphx_glr_auto_examples_image_compute_process_hires.py`
"""

import squidpy as sq

import matplotlib.pyplot as plt

###############################################################################
# First, we load an H&E stained tissue image.
# Here, we only load a cropped dataset to speed things up.
# In general, :func:`squidpy.im.process_img` can also process very large images
# (see :ref:`sphx_glr_auto_examples_image_compute_process_hires.py`).
img = sq.datasets.visium_hne_image_crop()

###############################################################################
# Then, we convert the image to grayscale and plot the result.
# With the argument ``img_id`` we can select the image layer that should be processed.
# When converting to grayscale, the channel dimensions change from 3 to 1.
# By default, the name of the resulting channel dimension will be ``'{{original_channel_name}}_gray'``.
# Use the argument ``channel_id`` to set a new channel name explicitely.
# By default, the resulting image is saved in the layer ``image_gray``.
# This behavior can be changed with the arguments ``copy`` and ``key_added``.

sq.im.process_img(img, img_id="image", processing="gray")

fig, axes = plt.subplots(1, 2)
axes[0].imshow(img["image"])
axes[0].set_title("original")
axes[1].imshow(img["image_gray"].squeeze(), cmap="gray")
axes[1].set_title("grayscale")
for ax in axes:
    ax.axis("off")
