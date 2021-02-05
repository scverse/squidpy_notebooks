"""
Use the image container
-----------------------

This example shows how to use :class:`squidpy.im.ImageContainer` to interact with image structured data.

"""

import squidpy as sq

import matplotlib.pyplot as plt

###############################################################################
# Load a pre-loaded image.
# To load your own data, use the image container constructor: 
# squidpy.im.ImageContainer()
img = sq.datasets.visium_hne_image()

###############################################################################
# Representation image in container:
# The image(s) are in the .data attribute of the instance, which is an
# xarray.Dataset. Note that this is a Dataset so that this attribute can hold
# multuple image-structured layers.
print(img.data)

###############################################################################
# You can access specific image-structured arrays in the imgage using their
# names.
print(img.data["image"])

###############################################################################
# Lazy loading:
# The image data can be lazily loaded why netcdf and explicitly loaded into
# memory via .load() and save to disk via .save():
img.load()

###############################################################################
# You can add images into the Data set via .add_img():
# Here we are adding the same image again under a different name as a toy example.
# It shares the same channel dimension with "image", so we can use the same
# label for "channel_id" here.
img.add_img(
    img=img.data["image"],
    img_id="image2",
    channel_id="channels",
    lazy=True,
)
