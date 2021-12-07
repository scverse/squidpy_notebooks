#!/usr/bin/env python
"""
Generate cropped images using :meth:`squidpy.im.ImageContainer.generate_spot_crops`
--------------------
This example shows how to use :meth:`squidpy.im.ImageContainer.generate_spot_crops`.

High-resolution tissue slides might be too large to fit in the memory.
Therefore, we use a generator that produces cropped images from the original image container object.
:meth:`squidpy.im.ImageContainer.generate_spot_crops` iterates over :attr:`anndata.AnnData.obsm` and extracts crops.
Implemented for 10X spatial datasets.
For Z-stacks, the specified ``library_id`` or list of ``library_id`` need to match the name of the Z-dimension.
Always extracts 2D crops from the specified Z-dimension.

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
# The arguments supported by :meth:`squidpy.im.ImageContainer.generate_spot_crops` are the following:
# * :class:`anndata.AnnData`, an annotated data object.
# * spatial_key (``str``), a key in ``anndata.AnnData.obsm`` where spatial coordinates are stored.
# * ``spot_scale``, a float, is scaling factor for the spot diameter. Larger values mean more context.
# * library_id (``Optional``[``str``]) –
#       *If _None_, there should only exist one entry in ``anndata.AnnData.uns`` [``'{spatial_key}'``].
#       *If a ``str``, first search ``anndata.AnnData.obs`` [``'{library_id}'``] which contains the mapping
#        from observations to library ids, then search ``anndata.AnnData.uns`` [``'{spatial_key}'``].
# * ``obs_names`` (``Optional``[``Iterable``[``Any``]]) – Observations from ``anndata.AnnData.obs_names``
#   for which to generate the crops. If _None_, all observations are used.
# * ``as_array``, if a ``str``, yields numpy array of the specified layer.
#   For example, ```as_array="image"``` specifies that we will crop from the "image" layer.
#   We can specify multiple layers to obtain crops from multiple pre-processing steps.
#       * ``as_array``, if ``True``, returns a ``dict`` object. Else, returns :class:`squidpy.im.ImageContainer`.
# * ``squeeze``, a ``bool`` removes singleton dimensions from the results, if ```as_array = True```.
# * ``return_obs``, a ``bool`` to yield names from ``obs_names``.
# * kwargs (``Any``), keyword arguments for ``crop_center()``.
# The ``gen`` can be looped on to generate the cropped images.
gen = img.generate_spot_crops(adata, scale=0.5, as_array="image", squeeze=True)

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
# We will now see how the cropped images differ with change in ``spot_size``.
# Scale=1 would crop the spot with exact coordinates.
# To illustrate, we change the spot_size and plot the images again by looping on ``next(gen)``.
gen = img.generate_spot_crops(adata, scale=1.5, as_array="image", squeeze=True)
fig, axes = plt.subplots(1, 5)
fig.set_size_inches((20, 6))
for i in range(5):
    axes[i].set_title(f"Cropped image {i+1}")
    axes[i].axis("off")
    axes[i].imshow(next(gen))
###############################################################################
# We can see the increase in the context with increase in the ``spot_size``.
gen = img.generate_spot_crops(adata, spot_scale=2, as_array="image", squeeze=True)
fig, axes = plt.subplots(1, 5)
fig.set_size_inches((20, 6))
for i in range(5):
    axes[i].set_title(f"Cropped image {i+1}")
    axes[i].axis("off")
    axes[i].imshow(next(gen))

###############################################################################
# Argument ``as_array`` also takes boolean ``True`` to return a ``dict`` where the keys are layers and
# values are ```numpy.ndarray```.
gen = img.generate_spot_crops(adata, spot_scale=0.5, as_array=True, squeeze=True)
for _i in range(1):
    print(next(gen))

###############################################################################
# Passing ``False`` to the argument ``as_array`` returns a :class:`squidpy.im.ImageContainer`.
gen = img.generate_spot_crops(adata, spot_scale=1.5, as_array=False, squeeze=True)
for _i in range(5):
    print(next(gen))

###############################################################################
# If ```return_obs = True```, yields a ``tuple`` (cropped image, ``obs_name``). Otherwise, yields just the crops.
# The type of the crops depends on ```as_array``` and the number of dimensions on ``squeeze``.
# This can be used to train the classifiers in machine learning.
gen = img.generate_spot_crops(adata, spot_scale=2, as_array=False, squeeze=True, return_obs=True)
for _i in range(5):
    print(next(gen))
