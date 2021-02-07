#!/usr/bin/env python
"""
Smoothing an image
------------------
This example shows how to use :func:`squidpy.im.process_img` to smooth an image layer of
:class:`squidpy.im.ImageContainer`.

We use the argument ``processing = 'smooth'`` to smooth the image. This calls :func:`skimage.filters.gaussian`
in the background. Keyword arguments ``kwargs`` are passed to the wrapped function.
This allows us to set the width of the gaussian kernel :math:`sigma`, used for smoothing.

.. seealso::

    - :ref:`sphx_glr_auto_examples_image_compute_gray.py`.
    - :ref:`sphx_glr_auto_examples_image_compute_process_hires.py`.
"""

import squidpy as sq

import matplotlib.pyplot as plt

# load H&E stained tissue image
img = sq.datasets.visium_hne_image_crop()

###############################################################################
# Smooth the image with ``sigma = 2``.
# With the argument ``img_id`` we can select the image layer that should be processed.
# By default, the resulting image is saved in the layer ``image_smooth``.
# This behavior can be changed with the arguments ``copy`` and ``key_added``.

sq.im.process_img(img, img_id="image", processing="smooth", sigma=2)

###############################################################################
# Now we can look at the result on a cropped part of the image.
crop = img.crop_corner(0, 0, 200, 200)

fig, axes = plt.subplots(1, 2)
axes[0].imshow(crop["image"])
axes[0].set_title("original")
axes[1].imshow(crop["image_smooth"])
axes[1].set_title("smoothed")
for ax in axes:
    ax.axis("off")
