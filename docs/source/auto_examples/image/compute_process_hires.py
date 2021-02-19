#!/usr/bin/env python
"""
Process a high-resolution image
-------------------------------

This example shows how to use :func:`squidpy.im.process` with tiling.

The function can be applied to any method (e.g., smoothing, conversion to grayscale)
or ``layer`` of a high-resolution image layer of :class:`squidpy.im.ImageContainer`.

By default, :func:`squidpy.im.process` processes the entire input image at once.
In the case of high-resolution tissue slides however, the images might be too big to fit in memory
and cannot be processed at once.
In that case you can use the argument ``size`` to tile the image in crops of shape ``size``,
process each crop, and re-assemble the resulting image.
Note that you can also use :func:`squidpy.im.segment` in this manner.

Note that depending on the processing function used, there might be border effects occurring at the edges
of the crops. In a future version, we will support the extraction of overlapping crops,
which can mitigate these effects.

.. seealso::

    - :ref:`sphx_glr_auto_examples_image_compute_smooth.py`.
    - :ref:`sphx_glr_auto_examples_image_compute_gray.py`.
    - :ref:`sphx_glr_auto_examples_image_compute_segment_fluo.py`.
"""

import squidpy as sq

import matplotlib.pyplot as plt

# load H&E stained tissue image
img = sq.datasets.visium_hne_image()

###############################################################################
# We will process the image by tiling it in crops of shape ``size = (1000, 1000)``.

sq.im.process(img, layer="image", method="gray", size=1000)

###############################################################################
# Now we can look at the result on a cropped part of the image.
crop = img.crop_corner(4000, 4000, size=2000)

fig, axes = plt.subplots(1, 2)
crop.show("image", ax=axes[0])
_ = axes[0].set_title("original")
crop.show("image_gray", cmap="gray", ax=axes[1])
_ = axes[1].set_title("grayscale")
