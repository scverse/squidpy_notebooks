#!/usr/bin/env python
r"""
Smooth an image
---------------

This example shows how to use :func:`squidpy.im.process` to smooth an image layer of :class:`squidpy.im.ImageContainer`.

We use the argument ``method="smooth"`` to smooth the image.
This calls :func:`skimage.filters.gaussian` in the background.
Keyword arguments ``kwargs`` are passed to the wrapped function.
This allows us to set the width of the Gaussian kernel, :math:`\\sigma`, used for smoothing.

.. seealso::

    - :ref:`sphx_glr_auto_examples_image_compute_gray.py`.
    - :ref:`sphx_glr_auto_examples_image_compute_process_hires.py`.
"""

import squidpy as sq

import matplotlib.pyplot as plt

# load the H&E stained tissue image
img = sq.datasets.visium_hne_image_crop()

###############################################################################
# Smooth the image with ``sigma = 2``.
# With the argument ``layer`` we can select the image layer that should be processed.
# By default, the resulting image is saved in the layer ``image_smooth``.
# This behavior can be changed with the arguments ``copy`` and ``layer_added``.
sq.im.process(img, layer="image", method="smooth", sigma=2)

###############################################################################
# Now we can look at the result on a cropped part of the image.
crop = img.crop_corner(0, 0, size=200)

fig, axes = plt.subplots(1, 2)
for i, layer in enumerate(["image", "image_smooth"]):
    crop.show(layer, ax=axes[i])
    axes[i].set_title(layer)
