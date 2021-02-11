#!/usr/bin/env python
"""
Crop images with ImageContainer
-------------------------------

This example shows how crop images from :class:`squidpy.im.ImageContainer`.

Specifically, it shows how to use:

    - :meth:`squidpy.im.ImageContainer.crop_corner()`
    - :meth:`squidpy.im.ImageContainer.crop_center()`

.. seealso::

    See :ref:`sphx_glr_auto_examples_image_compute_image_container.py` for general usage of
    :class:`squidpy.im.ImageContainer`.
"""

import squidpy as sq

import matplotlib.pyplot as plt

###############################################################################
# Load a H&E Visium image.
img = sq.datasets.visium_hne_image_crop()

###############################################################################
# Extracting single crops:
# Crops need to be sized and located. We distinguish crops located based on a
# corner coordinate of the crop and crops located based on the center coordinate
# of the crop.
# You can specify the crop coordinates in pixels (as ``int``) or in percentage of total image size (as ``float``)

crop_corner = img.crop_corner(1000, 1000, size=400)

crop_center = img.crop_center(1200, 1200, radius=200)

fig, axes = plt.subplots(1, 2)
crop_corner.show(ax=axes[0])
crop_center.show(ax=axes[1])

###############################################################################
# The result of the cropping functions is another ImageContainer

crop_corner
