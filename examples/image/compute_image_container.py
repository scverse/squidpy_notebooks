"""
ImageContainer object
---------------------

This example shows how to use :class:`squidpy.im.ImageContainer` to interact with image structured data.

.. seealso::

    See :ref:`sphx_glr_auto_examples_image_compute_crops.py` for examples how to crop images using
    :class:`squidpy.im.ImageContainer`.

"""

import squidpy as sq

###############################################################################
# Load a pre-loaded image.
# To load your own data, use the ImageContainer constructor:
# ``squidpy.im.ImageContainer(<image-path-or-array>)``
img = sq.datasets.visium_hne_image()

###############################################################################
# Representation image in container:
# The image(s) are in the :attr:`img.data` attribute of the instance, which is an
# :class:`xarray.Dataset`. Note that this is a Dataset so that this attribute can hold
# multiple image-structured layers.
print(img.data)

###############################################################################
# You can access specific image-structured arrays in the image using their
# names.
print(img["image"])

###############################################################################
# Lazy loading:
# The image data can be lazily loaded with `netcdf` and explicitly loaded into
# memory via ``.data.load()`` and saved to disk via ``.save()``:
img.data.load()

###############################################################################
# You can add images into the ImageContainer using ``.add_img()``:
# Here we are adding the same image again under a different name as a toy example.
# It shares the same channel dimension with "image", so we can use the same
# label for ``channel_dim`` here.
# If the added image layer has a different channel dimension, just specify a new
# label for ``channel_dim``.
img.add_img(
    img=img.data["image"],
    layer="image2",
    channel_dim="channels",
    lazy=True,
)
img
