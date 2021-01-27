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
Brief overview of Squidpy API
=============================

Import packages and load dataset
"""

import scanpy as sc
import squidpy as sq

adata = sq.datasets.visium_hne_adata()
img = sq.datasets.visium_hne_image()

###############################################################################
# Image analysis and feature extraction

# cell segmentation
sq.im.segment_img(img, img_id="image", model_group="watershed", thresh=70, geq=False)
# feature extraction
sq.im.calculate_image_features(
    adata,
    img,
    features=["summary", "segmentation"],
    features_kwargs={"segmentation": {"label_img_id": "segmented_watershed"}},
)

# plot mean intensity per spot (summary feature) and number of cells per spot (segmentation feature)
sc.pl.spatial(sq.pl.extract(adata, "img_features"), color=["summary_quantile_0.5_ch_0", "segmentation_label"])

###############################################################################
# Spatial graph analysis

# calculate graph
sq.gr.spatial_neighbors(adata)

# neighborhood enrichement
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster", figsize=(4, 3))

# spatially variable genes
sq.gr.moran(adata, n_jobs=6, n_perms=100)
sc.pl.spatial(adata, color=adata.uns["moranI"].index[:2])

# ligand-receptor interactions using cellphoneDB
sq.gr.ligrec(adata, cluster_key="cluster")
sq.pl.ligrec(
    adata,
    cluster_key="cluster",
    source_groups="Hippocampus",
    target_groups=["Pyramidal_layer", "Pyramidal_layer_dentate_gyrus"],
    means_range=(3, float("inf")),
    swap_axes=True,
)
