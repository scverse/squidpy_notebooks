# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: sphinx
#       format_version: '1.1'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3.9.12 ('squidpy_3_9')
#     language: python
#     name: python3
# ---

"""
# Feature extraction using CellProfiler

In this tutorial, we show how to use Squidpy with functions from [CellProfiler](https://cellprofiler.org) pipelines for image processing and feature extraction.   

We'll go through the following steps:

1. Load Visium fluorescence data.
2. Segment cells in Squidpy.
3. Calculate CellProfiler's granularity features for image crops of Visium spots.
4. Compute clustering on the image features in Squidpy.

CellProfiler is typically used via its GUI interface to build image processing pipelines.
First, download and install CellProfiler from the [download page](https://github.com/CellProfiler/CellProfiler). Check the issues on [CellProfiler Github](https://github.com/CellProfiler/CellProfiler/issues) in case of installation problems (can be tricky).   

Note: In the future, CellProfiler functions will also be accessible via Python directly (according to the announcements of the CellProfiler team). Since this is not yet publicly well documented, we'll restrict this tutorial to the laborious way of saving intermediate files to bridge Squidpy and CellProfiler. For information on how to use the `cellprofiler-core` package for Python integration, the following [link](https://github.com/CellProfiler/CellProfiler/wiki/CellProfiler-as-a-Python-package) might be helpful.   
"""

###############################################################################
# ## Import packages & data

import pandas as pd
import scanpy as sc
import anndata as ad
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from PIL import Image
import itertools

sc.logging.print_header()
print(f"squidpy=={sq.__version__}")
sc.set_figure_params(facecolor="white", figsize=(8, 8))

# load the pre-processed dataset
img = sq.datasets.visium_fluo_image_crop()
adata = sq.datasets.visium_fluo_adata_crop()

""
sq.pl.spatial_scatter(adata, color="cluster")

###############################################################################
# ## Setup files
# Create some arbitrary empty base directory `BASE_DIR` and a folder `imgs_dir = BASE_DIR + "tmp_imgs"` in it.

BASE_DIR = './squidpy_tutorial_cellprofiler/'
imgs_dir = BASE_DIR + 'tmp_imgs/'
Path(imgs_dir).mkdir(parents=True, exist_ok=True)

###############################################################################
# ## Image segmentation

sq.im.process(
    img=img,
    layer="image",
    method="smooth",
    sigma=[2,2,0,0],
)

sq.im.segment(img=img, layer="image_smooth", method="watershed", channel_ids=0, chunks=1000)

# plot the resulting segmentation
fig, ax = plt.subplots(1, 2)
img_crop = img.crop_corner(2000, 2000, size=500)
img_crop.show(layer="image", channel=0, ax=ax[0])
img_crop.show(
    layer="segmented_watershed",
    channel=0,
    ax=ax[1],
)

###############################################################################
# ## Save image crops and segmentations
# For each Visium spot we generate crops:

for crop,obs in itertools.islice(img.generate_spot_crops(adata, obs_names=adata.obs_names,
                                                         return_obs=True, as_array=False), 2):
    fig, axes = plt.subplots(1, 2)
    crop.show('image_smooth', ax=axes[0])
    crop.show('segmented_watershed', ax=axes[1], cmap='Greys')
    plt.show() 

###############################################################################
# Save the fluorescence image and the segmentations of each crop.

for crop, obs in img.generate_spot_crops(adata, obs_names=adata.obs_names,
                                         return_obs=True, as_array=False):
    Image.fromarray((np.array(crop['image_smooth']) * 255).squeeze().astype(np.uint8)).save(imgs_dir + f'image_{obs}.png')
    Image.fromarray((np.array(crop['segmented_watershed'][:, :, 0])).squeeze()).save(imgs_dir + f'segments_{obs}.tif')

###############################################################################
# ## CellProfiler Pipeline: Calculate Image Features
#
# ### 1. Open the CellProfiler GUI App. A new project should open automatically.
# ![image1](tutorial_cellprofiler_images/01CP_new_project.png)
#
# ### 2. Save the project in `BASE_DIR`.
# ![image2](tutorial_cellprofiler_images/02CP_save_project.png)
#
# ### 3. CP-Pipeline `Images`: Drag and Drop the folder `imgs_dir` into CellProfiler.
# ![image3](tutorial_cellprofiler_images/03CP_import_images.png)
#
# ### 4. CP-Pipeline `NamesAndTypes`: Declare to load crops as color images and segmentations as objects. Crops and segmentations files are aligned automatically.
# ![image4](tutorial_cellprofiler_images/04CP_NamesAndTypes.png)
#
# ### 5. Convert image crops to gray images: Add `ColorToGray` module and define parameters.
# ![image5](tutorial_cellprofiler_images/05CP_Create_ColorToGray.png)
# ![image6](tutorial_cellprofiler_images/06CP_ColorToGray.png)
#
# ### 6. Measure CellProfiler's Granularity features within segments for each crop: Add `MeasureGranularity` module and define parameters.
# ![image7](tutorial_cellprofiler_images/07CP_Create_Granularity.png)
# ![image8](tutorial_cellprofiler_images/08CP_Granularity.png)
#
# ### 7. Export results to csv.
# ![image9](tutorial_cellprofiler_images/09CP_Create_Export.png)   
# ![image10](tutorial_cellprofiler_images/10CP_Export.png)
#
# ### 8. Run CellProfiler Pipeline (takes several minutes).
# ![image11](tutorial_cellprofiler_images/11CP_Run.png)

###############################################################################
# You can now delete the images in `imgs_dir` by uncommenting the cell below:

# # !rm -rf $imgs_dir

###############################################################################
# ## Cluster features
# Load CellProfiler output of Granularity features of each Visium spot and rearrange the output for `adata.obsm`.

# load CellProfiler output
df = pd.read_csv(BASE_DIR + '/squidpy_Image.csv')
# set obs names as index
df.index = df['FileName_Spot_Crop'].apply(lambda s: s.split("_")[1].split(".")[0])
df.index.name = "obs"
# Get the measured Granularity features and rename them
features = [f'{stat}_Segments_Granularity_{i}_SpotGray' for stat in ['Mean', 'StDev'] for i in range(1, 17)]
df = df[features]
df.columns = [s.split("_")[0] + s.split("_")[2] + s.split("_")[3] for s in df.columns]

""
df.head()


###############################################################################
# Load image features into `adata.obsm["features"]` and compute Leiden clustering.

def cluster_features(features: pd.DataFrame):
    """Calculate leiden clustering of features.
    """
    # create temporary adata to calculate the clustering
    adata = ad.AnnData(features)
    # important - feature values are not scaled, so need to scale them before PCA
    sc.pp.scale(adata)
    # calculate leiden clustering
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)

    return adata.obs["leiden"]

# Save feature in adata.obsm
adata.obsm["features"] = df.loc[adata.obs_names].copy()

# Calculate leiden clustering
adata.obs["features_granularity_cluster"] = cluster_features(adata.obsm["features"])

# Plot clusterings
sq.pl.spatial_scatter(
    adata,
    color=[
        "features_granularity_cluster",
        "cluster",
    ],
    ncols=2,
)

###############################################################################
# We see that e.g. the blue cluster is prominent in the Hippocampus and the Thalamus, and the brown cluster is prominent in dense regions like Dentate gyrus.   
#
# This tutorial is an example on how to integrate CellProfiler within Squidpy pipelines. You can adapt this scheme to other task distributions: e.g. segmentation with CellProfiler and image feature calculation in Squidpy, or process images with CellProfiler, followed by a Squidpy pipeline etc.
