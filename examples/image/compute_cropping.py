"""
Use the image container for cropping images
-------------------------------------------

This example shows how to use:

    - :func:`squidpy.im.ImageContainer.crop_corner()`
    - :func:`squidpy.im.ImageContainer.crop_center()`

"""

import squidpy as sq

import matplotlib.pyplot as plt

###############################################################################
# Load a pre-loaded image.
img = sq.datasets.visium_hne_image()

###############################################################################
# Extracting single crops:
# Crops need to be sized and located. We distinguish crops located based on a
# corner coordinate of the crop and crops located based on the centre coordinate
# of the crop:

crop = img.crop_corner(0, 0, 200, 200)

fig, axes = plt.subplots(1, 2)
axes[0].imshow(crop["image"])
axes[0].set_title("original")
axes[1].imshow(crop["image_smooth"])
axes[1].set_title("smoothed")
for ax in axes:
    ax.axis("off")

crop = img.crop_center(100, 100, 200, 200)

fig, axes = plt.subplots(1, 2)
axes[0].imshow(crop["image"])
axes[0].set_title("original")
axes[1].imshow(crop["image_smooth"])
axes[1].set_title("smoothed")
for ax in axes:
    ax.axis("off")
