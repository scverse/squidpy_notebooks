# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: sphinx
#       format_version: '1.1'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

"""
# Nuclei segmentation using StarDist

In this tutorial, we show how we can use the `StarDist` segmentation method in `squidpy.im.segment` for nuclei segmentation. 

**StarDist** <cite data-cite="stardist_schmidt2018">Schmidt et al. (2018)</cite> and <cite data-cite="stardist_weigert2020">Weigert et al. (2020)</cite> , ([code](https://github.com/stardist/stardist))  uses star-convex polygons to localize cell for which a convolutional neural network was trained to predict pixel-wise polygons for each cell position. 

To run the notebook locally, create a conda environment as *conda env create -f stardist_environment.yml* using this [stardist_environment.yml](https://github.com/scverse/squidpy_notebooks/blob/main/envs/stardist_environment.yml), which installs Squidpy, TensorFlow, and StarDist.

**Note:** We frequently recognized a dying notebook kernel when importing other packages before StarDist with the following message "The kernel appears to have died. It will restart automatically." We therefore recommend to import StarDist first.
"""

# Import the StarDist 2D segmentation models.
from stardist.models import StarDist2D 
# Import the recommended normalization technique for stardist.
from csbdeep.utils import normalize

# Import squidpy and additional packages needed for this tutorial.
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# StarDist has four pre-trained models for 2D images. We will show an example for the Versatile (fluorescent nuclei) model and the Versatile (H&E nuclei). To use the StarDist model, we define a wrapper that normalizes the image with the recommended method, initializes the model and returns the segmentation mask. 

StarDist2D.from_pretrained()

###############################################################################
# The method parameter of the `sq.im.segment` method accepts any callable with the signature:
# `numpy.ndarray` ``(height, width, channels)`` **->** `numpy.ndarray` ``(height, width[, channels])``. Additional model specific arguments will also be passed on. 

###############################################################################
# ## Cell segmentation on Visium fluorescence data
#

# Load the image and visualize its channels.
img = sq.datasets.visium_fluo_image_crop()
crop = img.crop_corner(1000, 1000, size=1000)
crop.show(channelwise=True)

###############################################################################
# Additionally, we will have a look at the pre-trained StarDist model. The `2D_versatile_fluo model` works on one channel, as `n_channel_in = 1`. We will run the segmentation on the first channel of the image in this example. 

StarDist2D.from_pretrained('2D_versatile_fluo')


###############################################################################
# The input image is normalized beforehand by supplying a `normalizer` to the prediction function. We pass the recommended StarDist normalization method from `csbdeep.utils` into our callable.
#
# Calling `model.predict_instances` will:
#
# - predict object probabilities and star-convex polygon distances.
# - perform non-maximum suppression (with overlap threshold `nms_thresh`) for polygons above object probability threshold `prob_thresh`.
# - render all remaining polygon instances in a label image.
# - return the label instances image and also the details (coordinates, etc.) of all remaining polygons.
#
# For our purpose, we will only return the respective labels. Check the detailed example StarDist [notebook](https://nbviewer.jupyter.org/github/stardist/stardist/blob/master/examples/2D/3_prediction.ipynb) for more information.

def stardist_2D_versatile_fluo(img, nms_thresh=None, prob_thresh=None):
    # Make sure to normalize the input image beforehand or supply a normalizer to the prediction function.
    # this is the default normalizer noted in StarDist examples.
    img = normalize(img, 1, 99.8, axis=(0,1))
    model = StarDist2D.from_pretrained('2D_versatile_fluo')
    labels, _ = model.predict_instances(img, nms_thresh=nms_thresh, prob_thresh=prob_thresh)
    return labels


""
sq.im.segment(
    img=crop, 
    layer="image", 
    channel=0, 
    method=stardist_2D_versatile_fluo, 
    layer_added='segmented_stardist',
    nms_thresh=None,
    prob_thresh=None
)

""
# Plot the DAPI channel of the image crop and the segmentation result.
print(crop)
print(f"Number of segments in crop: {len(np.unique(crop['segmented_stardist']))}")

fig, axes = plt.subplots(1, 2)
crop.show("image", channel=0, ax=axes[0])
_ = axes[0].set_title("DAPI")
crop.show("segmented_stardist", cmap="jet", interpolation="none", ax=axes[1])
_ = axes[1].set_title("segmentation")

###############################################################################
# ## Cell segmentation on H&E stained tissue data

# load H&E stained tissue image and crop to a smaller segment
img = sq.datasets.visium_hne_image_crop()
crop = img.crop_corner(0, 0, size=1000)
crop.show("image")


""
def stardist_2D_versatile_he(img, nms_thresh=None, prob_thresh=None):
    #axis_norm = (0,1)   # normalize channels independently
    axis_norm = (0,1,2) # normalize channels jointly
    # Make sure to normalize the input image beforehand or supply a normalizer to the prediction function.
    # this is the default normalizer noted in StarDist examples.
    img = normalize(img, 1, 99.8, axis=axis_norm)
    model = StarDist2D.from_pretrained('2D_versatile_he')
    labels, _ = model.predict_instances(img, nms_thresh=nms_thresh, prob_thresh=prob_thresh)
    return labels


###############################################################################
# StarDist H&E segmentation method works on three input channels as `n_channel_in = 3`. We therefore pass `channel = None` 
#  to the `sq.img.segment` method which will then run the given segmentation method on all given channels.

StarDist2D.from_pretrained('2D_versatile_he')

""
sq.im.segment(
    img=crop, 
    layer="image", 
    channel=None,
    method=stardist_2D_versatile_he, 
    layer_added='segmented_stardist_default', 
    prob_thresh=None,
    nms_thresh=None
)

""
print(crop)
print(f"Number of segments in crop: {len(np.unique(crop['segmented_stardist_default']))}")

fig, axes = plt.subplots(1, 2)
crop.show("image", ax=axes[0])
_ = axes[0].set_title("H&H")
crop.show("segmented_stardist_default", cmap="jet", interpolation="none", ax=axes[1])
_ = axes[1].set_title("segmentation")

###############################################################################
# Adjusting the `prob_thresh` parameter will enhance the segmentation. We show this additionally for `prob_thresh = 0.3`. Please be aware that the print statement of the default values will remain unchanged, even if you adjusted the parameters.

sq.im.segment(
    img=crop, 
    layer="image", 
    channel=None,
    method=stardist_2D_versatile_he, 
    layer_added='segmented_stardist', 
    prob_thresh=0.3,
    nms_thresh=None
)

""
print(crop)
print(f"Number of segments in crop: {len(np.unique(crop['segmented_stardist']))}")

fig, axes = plt.subplots(1, 2)
crop.show("image", ax=axes[0])
_ = axes[0].set_title("H&H")
crop.show("segmented_stardist", cmap="jet", interpolation="none", ax=axes[1])
_ = axes[1].set_title("segmentation")
