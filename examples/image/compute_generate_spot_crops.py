#!/usr/bin/env python
"""
Generate cropped images using :meth:`squidpy.im.ImageContainer.generate_spot_crops`
--------------------
This example shows how to use :meth:`squidpy.im.ImageContainer.generate_spot_crops`.

High-resolution tissue slides might be too large to fit in the memory.
Therefore, we use a generator that produces cropped images from the original image container object.
:meth:`squidpy.im.ImageContainer.generate_spot_crops` iterates over :attr:`anndata.AnnData.obsm` and extracts crops.

.. seealso::
    - :ref:`sphx_glr_auto_examples_image_compute_crops.py`
    - :ref:`sphx_glr_auto_examples_image_compute_process_hires.py`
    - :ref:`sphx_glr_auto_examples_image_compute_gray.py`
"""

import squidpy as sq

import matplotlib.pyplot as plt

###############################################################################
# First, we load the H&E stained tissue image.
# Here, we only load a cropped dataset to speed things up.
# In general, :meth:`squidpy.im.ImageContainer.generate_spot_crops` can also process very large images.
# See :ref:`sphx_glr_auto_examples_image_compute_process_hires.py`.
# Second, we load the related anndata for the H&E stained tissue image.
img = sq.datasets.visium_hne_image_crop()
adata = sq.datasets.visium_hne_adata_crop()

###############################################################################
# Next, we use :meth:`squidpy.im.ImageContainer.generate_spot_crops` to make a generator that generates cropped images.
# In addition to :class:`anndata.AnnData` as a parameter, we can specify ``spot_scale``,``as_array``,``squeeze``.
# * :class:`anndata.AnnData` is an annotated data object.
# * ``spot_scale``, a float, is scaling factor for the spot diameter. Larger values mean more context.
# * ``as_array``, if string, yields numpy array of the specified layer.
# * ``squeeze``, if true, removes singleton dimensions from the results.
# The ``gen`` can be looped on to generate the cropped images.
gen = img.generate_spot_crops(adata, scale=2, as_array="image", squeeze=True)

###############################################################################
# When called, the ``next(gen)`` produces consecutive cropped images each time.
# Now, we plot the cropped images using matplotlib.
fig, axes = plt.subplots(1, 5)
fig.set_size_inches((20, 6))
for i in range(5):
    axes[i].set_title(f"Cropped image {i+1}")
    axes[i].axis("off")
    axes[i].imshow(next(gen))

###############################################################################
# To illustrate successive cropped images, we plot the images again by looping on ``next(gen)``.
fig, axes = plt.subplots(1, 5)
fig.set_size_inches((20, 6))
for i in range(5):
    axes[i].set_title(f"Cropped image {i+1}")
    axes[i].axis("off")
    axes[i].imshow(next(gen))
