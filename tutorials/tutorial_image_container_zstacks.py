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
In this example we showcase how to use z-stacks with :class:`squidpy.im.ImageContainer`

It is possible to acquire several consecutive image slices from the same tissue.
Squidpy's `ImageContainer` supports storing, processing, and visualization of these z-stacks.

Here, we use the Visisum 10x mouse brain sagittal slices as an example of a z-stack image with two Z dimensions.
We will use the "hires" images contained in the :class:`anndata.AnnData` object, but you could also use the
original resolution tiff images in the `ImageContainer`.

.. seealso::

    See :ref:`sphx_glr_auto_tutorials_tutorial_image_container.py` for a general introduction to the `ImageContainer`.


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
# To allow mapping from observations in `adata` to the correct Z dimension in `img`,
# we will store a ``library_id`` column in ``adata.obs`` and associate each ``library_id``
# to a Z dimension in the `ImageContainer`.
#
# For this, we will use :func:`anndata.concat` with ``uns_merge=only``
# (to ensure that `uns` entries are correctly concatenated),
# ``label='library_id'`` and ``keys=library_ids`` (to create the necessary column in ``adata.obs``.
#
# To concatenate the individual `ImageContainers`, we will use :meth:`squidpy.im.ImageContainer.concat`, specifying
# ``library_ids=library_ids`` for associating each image with the correct observations in `adata`.

adata = ad.concat(adatas, uns_merge="only", label="library_id", keys=library_ids, index_unique="-")
img = sq.im.ImageContainer.concat(imgs, library_ids=library_ids)

###############################################################################
# `adata` now contains a ``library_id`` column in ``adata.obs``, which maps observations to a unique `library_id`

print(adata)
adata.obs

###############################################################################
# `img` contains the 2D images concatenated along the Z dimension in one image layer.
# The Z dimensions are named the same as the `library_id`s in `adata` to allow a mapping from `adata` to `img`.

print(img["image"].z)
img

###############################################################################
# It is also possible to initialise the `ImageContainer` with images that already contain the Z dimension.
# In this case you need to specify the ``library_id`` argument in the constructor.
# In addition, you might want to set ``dims`` to the correct ordering of dimensions manually for more control.

arr = img["image"].values
print(arr.shape)
img2 = sq.im.ImageContainer(arr, library_id=library_ids, dims=("y", "x", "z", "channels"))
img2

###############################################################################
# Generally, an `ImageContainer` with more than one Z dimension can be used in the same way as an `ImageContainer`
# with only one Z dimension.
# In addition, we can specify `library_id` to cropping, preprocessing,
# and segmentation functions if we'd like to only process a specific `library_id`.

###############################################################################
# Visualization
# -------------
#
# For using `sc.pl.spatial`, subset the `adata` to the desired `library_id`

library_id = library_ids[0]
sc.pl.spatial(adata[adata.obs["library_id"] == library_id], library_id=library_id, color="in_tissue")

###############################################################################
# :meth:`squidpy.im.ImageContainer.show` works with z-stacks out of the box, by plotting them as separate images.
# Additionally, you can specify a `library_id` if you only want to plot one Z dimension.

img.show()

###############################################################################
# Interactive visualisation of z-stacks is also possible.
# The Napari viewer will have a slider at the bottom, allowing you to choose the Z dimension to display.
# The `adata` observations are automatically updated to the current Z dimension.
#
# When calling ``img.interactive`` just specify ``library_key`` as the column name in ``adata.obs``
# which maps from observations to `library_ids`
#
# .. code-block:: python
#
#     import napari
#     with napari.gui_qt():
#         img.interactive(adata, library_key='library_id')

###############################################################################
# Cropping
# --------
#
# By default, the cropping functions will crop all Z dimensions.

crop = img.crop_corner(500, 1000, size=500)
crop.show()


###############################################################################
# You can also specify ``library_id``, as either a single or multiple Z dimensions to crop.

img.crop_corner(500, 1000, size=500, library_id=library_ids[0]).show()

###############################################################################
# Processing and segmenting
# -------------------------
#
# Let us smooth the image.
# When not specifying a `library_id`, :func:`squidpy.im.process` treats the image as a 3D volume.
# As we would like to smooth only in x and y dimensions, and not in z, we need so specify a per-dimension `sigma`.
# The internal dimensions of the image are ``y,x,z,channels``, as you can check with ``crop['image'].dims``.
# Therefore, to only smooth in x and y, we need to specify ``sigma=[10,10,0,0]``.

sq.im.process(img, layer="image", method="smooth", sigma=[10, 10, 0, 0], layer_added="smooth1")
img.show("smooth1")

###############################################################################
# Now, let us just smooth one `library_id`.
# Specifying `library_id` means that the processing function will process each Z dimension separately.
# This means that now the dimensions of the processed image are ``y,x,channels`` (with ``z`` removed), meaning that
# we have to update `sigma` accordingly.
# If the number of channels does not change due to the processing, :func:`squidpy.im.process` implies the identity
# function for non-processed Z dimensions.

sq.im.process(img, layer="image", method="smooth", sigma=10, layer_added="smooth2", library_id=library_ids[0])
img.show("smooth2")

###############################################################################
# None, only the first `library_id` is smoothed.
# For the second, the original image was used.
#
# If the processing function changes the number of dimensions, non-processed Z dimensions will contain 0.
# Lets see this behavior with using ``method='gray'``, which moves from 3 channels (RGB) to one channel (gray).

sq.im.process(img, layer="image", method="gray", layer_added="gray", library_id=library_ids[0])
img.show("gray", cmap="gray")

###############################################################################
# :func:`squidpy.im.segment` works in the same way, just specify `library_id` if you only wish to
# segment specific Z dimensions.
#
# Feature calculation
# -------------------
#
# Calculating features from z-stack images is straight forward as well.
# With more than one Z dimension, we just need to specify the column name in ``adata.obs``
# which contains the mapping from observations to `library_ids`
# to allow the function to extract the features from the correct Z dimension.
# As of now, features can only be extracted on 2D, meaning from the Z dimension that the current spot is located on.
#
# The following call extracts features for each observation in `adata`, automatically choosing the correct
# Z dimension in `img`.

adata_crop = crop.subset(adata)  # subset adata to the image crop
sq.im.calculate_image_features(adata_crop, crop, library_id="library_id", layer="image", features="summary", n_jobs=4)
adata_crop.obsm["img_features"]


###############################################################################
# The calculated features can now be used in downstream Scanpy analyses, by e.g. using all Z dimensions
# to cluster spots based on image features and gene features.
#
# Here, we cluster genes and calculated features using a standard Scanpy workflow.

sc.pp.normalize_total(adata_crop, inplace=True)
sc.pp.log1p(adata_crop)
sc.pp.pca(adata_crop)
sc.pp.neighbors(adata_crop)
sc.tl.leiden(adata_crop)

sc.pp.neighbors(adata_crop, use_rep="img_features", key_added="neigh_features")
sc.tl.leiden(adata_crop, neighbors_key="neigh_features", key_added="leiden_features")

###############################################################################
# Visualize the result interactively using Napari, or statically using :func:`scanpy.pl.spatial`:
#
# .. code-block:: python
#
#     import napari
#     with napari.gui_qt():
#         img.interactive(adata, library_key='library_id')

sc.pl.spatial(
    adata_crop[adata_crop.obs["library_id"] == library_ids[0]],
    library_id=library_ids[0],
    color=["leiden", "leiden_features"],
)

sc.pl.spatial(
    adata_crop[adata_crop.obs["library_id"] == library_ids[1]],
    library_id=library_ids[1],
    color=["leiden", "leiden_features"],
)
