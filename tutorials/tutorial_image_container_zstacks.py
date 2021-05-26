# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: sphinx
#       format_version: '1.1'
#       jupytext_version: 1.9.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

"""
Use z-stacks with ImageContainer
================================
In this example we showcase how to use z-stacks with ImageContainer

It is possible to acquire several consecutive image slices from the same tissue.
Squidpy's ImageContainer supports storing, processing, and visualization of these z-stacks.

Here, we use the Visisum 10x mouse brain sagittal slices as an example of a z-stack image with two z dimensions.
We will use the "hires" images contained in `adata`, but you could also use the original resolution tiff
images in the ImageContainer.

Import Libraries and load individual image sections
---------------------------------------------------
"""

import scanpy as sc
import anndata as ad
import squidpy as sq

library_ids = ["V1_Mouse_Brain_Sagittal_Posterior", "V1_Mouse_Brain_Sagittal_Posterior_Section_2"]

adatas = []
imgs = []
use_hires_tiff = False
for library_id in library_ids:
    adatas.append(sc.datasets.visium_sge(library_id, include_hires_tiff=use_hires_tiff))
    adatas[-1].var_names_make_unique()
    if use_hires_tiff:
        imgs.append(sq.im.ImageContainer(adatas[-1].uns["spatial"][library_id]["metadata"]["source_image_path"]))
    else:
        # as we are using a scaled image, we need to specify a scalefactor
        # to allow correct mapping to adata.obsm['spatial']
        imgs.append(
            sq.im.ImageContainer(
                adatas[-1].uns["spatial"][library_id]["images"]["hires"],
                scale=adatas[-1].uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"],
            )
        )

###############################################################################
# Concatenate per-section data to a z-stack
# -----------------------------------------
#
# To allow mapping from observations in the `adata` to the correct z dimension in `img`,
# we will store a `library_id` column in `adata.obs` and associate each `library_id`
# to a z dimension in the `ImageContainer`.
#
# For this, we will use `ad.concat` with ``uns_merge=only`` (to ensure that `uns` entries are correctly concatenated),
# `label='library_id'` and `keys=library_ids` (to create the necessary `adata.obs`.
# To concatenate the individual ImageContainers, we just need to specify `library_ids=library_ids`
# to `sq.im.ImageContainer.concat`.

adata = ad.concat(adatas, uns_merge="only", label="library_id", keys=library_ids, index_unique="-")
img = sq.im.ImageContainer.concat(imgs, library_ids=library_ids)

###############################################################################
# adata now contains a `library_id` column in `adata.obs`, which maps observations to a unique `library_id`

print(adata)
adata.obs

###############################################################################
# img contains the 2D images concatenated along the z dimension.
# The z dimensions are named like the `library_id`s in `adata` to allow a mapping from `adata` to `img`.

print(img["image"].z)
img

###############################################################################
# Generally, an `ImageContainer` with more than one z dim can be used in the same way as an `ImageContainer`
# with only one z dimension.
# In addition, we can specify `library_id` to cropping, preprocessing,
# and segmentation functions if we'd like to only process a specific `library_id`.

###############################################################################
# Plot AnnData and ImageContainer
# -------------------------------
#
# For using `sc.pl.spatial`, subset the `adata` to the desired `library_id`

library_id = library_ids[0]
sc.pl.spatial(adata[adata.obs["library_id"] == library_id], library_id=library_id, color="in_tissue")

###############################################################################
# `img.show` works with z-stacks out of the box, by plotting them as separate images.
# Additionally, you can specify a `library_id` if you only want to plot one z dimension.

img.show()

###############################################################################
# Cropping z-stacks
# -----------------

crop = img.crop_corner(500, 1000, size=500)
crop.show()

###############################################################################
# when cropping, you can specify the `library_id`.
# This argument takes either a string (single library_id) or a list of library_ids to crop

img.crop_corner(500, 1000, size=500, library_id=library_ids[0])

###############################################################################
# Processing and segmenting z-stacks
# ----------------------------------
#
# let us smooth the image.
# when not specifying a `library_id`, the processing function will treat the image as a 3D volume.
# As we would like to smooth only in x and y dimensions, and not in z, we need so specify a per-dimension `sigma`.
# The internal dimensions of the image are y,x,z,channels, as you can check with `crop['image'].dims`.
# Therefore, to only smooth in x and y, we need to specify sigma=[10,10,0,0].

sq.im.process(img, layer="image", method="smooth", sigma=[10, 10, 0, 0], layer_added="smooth1")
img.show("smooth1")

###############################################################################
# now, let us just smooth one `library_id`.
# Specifying `library_id` means that the processing function will process each z dimension separately.
# This means that now the dimensions of the processed image are y,x,channels (with z removed), meaning that
# we have to update `sigma` accordingly.
# If the number of channels does not change due to the processing, `process` will use the original image
# for non-processed z dimensions

sq.im.process(img, layer="image", method="smooth", sigma=[10, 10, 0], layer_added="smooth2", library_id=library_ids[0])
img.show("smooth2")

###############################################################################
# only the first library_id is smoothed. For the second, the original image was placed.
#
# if the processing function changes the number of dimensions, non-processed z dimensions will contain 0.
# Lets see this behavior with using method='gray', which moves from 3 channels (RGB) to one channel (gray).

sq.im.process(img, layer="image", method="gray", layer_added="gray")
img.show("gray", cmap="gray")

###############################################################################
# `segment` works in the same way, just specify `library_id` if you only wish to segment specific z dimensions.
#
# Calculate features from z-stacks
# --------------------------------
# Calculating features from z-stack images is straight forward as well.
# With more than one z dimension, the supplied `adata` needs to contain a `library_id` column on `adata.obs`
# to allow the function to extract the features from the correct z dimension.
# As of now, features can only be extracted on 2D, meaning from the z dimension that the current spot is located on.
#
# The following call extracts features for each obs in `adata`, automatically choosing the correct z dimension in `img`.

adata_crop = crop.subset(adata)
sq.im.calculate_image_features(adata_crop, crop, layer="image", features="summary", n_jobs=4)
adata_crop.obsm["img_features"]


###############################################################################
# TODO: interactive visualisation of z-stacks
# -------------------------------------------
