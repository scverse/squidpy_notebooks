"""
Processing a high-resolution Image
----------------------------------

This example shows how to use :func:`squidpy.im.process_img` to apply any processing function
(smoothing, conversion to grayscale) to a high-resolution image layer of :class:`squidpy.im.ImageContainer`.

By default, :func:`squidpy.im.process_img` processes the entire input image at once.
In the case of high-resolution tissue slides however, the images might be too big to fit in memory
and cannot be processed at once.
In that case you can use the arguments ``xs`` and ``ys`` that will tile the image in crops of size ``(ys, xs)``,
process each crop, and re-assemble the resulting image.
Note that you can also use :func:`squidpy.im.segment_img` in this manner.

Note that depending on the processing function used, there might be border effects occurring at the edges
of the crops. In a future version, we will support the extraction of overlapping crops,
which can mitigate these effects.

For more usage examples see also   :ref:`sphx_glr_auto_examples_image_compute_smooth.py`,
:ref:`sphx_glr_auto_examples_image_compute_gray.py`, and
:ref:`sphx_glr_auto_examples_image_compute_segment_fluo.py`.
"""

import squidpy as sq

import matplotlib.pyplot as plt

# load H&E stained tissue image
img = sq.datasets.visium_hne_image()

###############################################################################
# We will process the image by tiling it in crops of shape ``(ys, xs) = (1000, 1000)``.

sq.im.process_img(img, img_id="image", processing="gray", xs=1000, ys=1000)

###############################################################################
# Now we can look at the result on a cropped part of the image.
crop = img.crop_corner(4000, 4000, 2000, 2000)

fig, axes = plt.subplots(1, 2)
axes[0].imshow(crop["image"])
axes[0].set_title("original")
axes[1].imshow(crop["image_gray"].squeeze(), cmap="gray")
axes[1].set_title("converted to grayscale")
for ax in axes:
    ax.axis("off")
