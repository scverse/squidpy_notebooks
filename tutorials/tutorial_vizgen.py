# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: sphinx
#       format_version: '1.1'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

"""
Analyze Vizgen data
===================
"""

import scanpy as sc
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path

sc.logging.print_header()

###############################################################################
# Download the data from [Vizgen MERFISH Mouse Brain Receptor Dataset](https://info.vizgen.com/mouse-brain-map?submissionGuid=a66ccb7f-87cf-4c55-83b9-5a2b6c0c12b9). Unpack the `.tar.gz` file. The dataset contains a MERFISH measurement of a gene panel containing 483 total genes including canonical brain cell type markers, GPCRs, and RTKs measured on 3 full coronal slices across 3 biological replicates. This is one slice of replicate 1.
#
# Unfortunately, the data needs to be downloaded manually. You need these 3 files in a new folder `tutorial_data` in the same path as your notebook. 
# - `datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_by_gene_S1R1.csv`
# - `datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_metadata_S1R1.csv`
# - `datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_images_micron_to_mosaic_pixel_transform.csv`
# The last file should be in the `images` folder.
#
# The following lines create the folder structure which can be use to load the data.

# # # Download and unpack the Vizgen data
# # !mkdir tutorial_data
# # !mkdir tutorial_data/vizgen_data
# # !mkdir tutorial_data/vizgen_data/images

""
vizgen_dir = Path().resolve() / "tutorial_data" / "vizgen_data"

adata = sq.read.vizgen(
    path=vizgen_dir,
    counts_file="datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_by_gene_S1R1.csv",
    meta_file="datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_metadata_S1R1.csv",
    transformation_file="datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_images_micron_to_mosaic_pixel_transform.csv",
)

###############################################################################
# Make the variable names unique using the method `anndata.var_names_make_unique`.
# Obtain the mitochondrial genes using their names prefixed with "mt-".
# Calculate the quality control metrics on the `anndata.AnnData` using `scanpy.pp.calculate_qc_metrics`.

adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"],percent_top=(50, 100, 200, 300), inplace=True)

###############################################################################
# Calculate quality control metrics
# ---------------------------------

###############################################################################
# Plot the "total_counts" and "n_genes_by_counts" from `anndata.obs`.
# The first subplot shows `adata.obs["total_counts"]`, the second `adata.obs["total_counts"]` less than 10000.
# The third subplot displays `adata.obs["n_genes_by_counts"]` while the fourth displays the `adata.obs["n_genes_by_counts"]` less than 4000.

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.distplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)
sns.distplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1],
)
sns.distplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    bins=60,
    ax=axs[2],
)
sns.distplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3],
)

###############################################################################
# Filter the cells based on the minimum number of counts required using `scanpy.pp.filter_cells`. Filter the genes based on the minimum number of cells required with `scanpy.pp.filter_genes`. The parameters for the both were specified based on the plots above. They were set to filter out the cells and genes with minimum counts and minimum cells respectively.

sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_genes(adata, min_cells=10)

###############################################################################
# Annotate the highly variable genes based on the count data by using `scanpy.pp.highly_variable_genes` with `flavor="seurat_v3"`. Normalize counts per cell using `scanpy.pp.normalize_total`.
#
# Logarithmize, do principal component analysis, compute a neighborhood graph of the observations using `scanpy.pp.log1p`, `scanpy.pp.pca` and `scanpy.pp.neighbors` respectively.
#
# Use `scanpy.tl.umap` to embed the neighborhood graph of the data and cluster the cells into subgroups employing `scanpy.tl.leiden`.

###############################################################################
# You may have to install `scikit-misc` package for highly variable genes identification.

# # !pip install scikit-misc

""
adata.layers["counts"] = adata.X.copy()
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

###############################################################################
# Visualize annotation on UMAP and spatial coordinates
# ----------------------------------------------------

###############################################################################
# Subplot with scatter plot in UMAP (Uniform Manifold Approximation and Projection) basis. The embedded points were colored, respectively, according to the total counts, number of genes by counts, and leiden clusters in each of the subplots. This gives us some idea of what the data looks like.

sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,
)

""
sq.pl.spatial_scatter(adata, shape=None,color=[
        "leiden",
    ],
    wspace=0.4,)

###############################################################################
# ----

###############################################################################
# Computation of spatial statistics
# ---------------------------------

###############################################################################
# Building the spatial neighbors graphs
# -------------------------------------

###############################################################################
# This example shows how to compute centrality scores, given a spatial graph and cell type annotation.
#
# The scores calculated are closeness centrality, degree centrality and clustering coefficient with the following properties:
# * closeness centrality - measure of how close the group is to other nodes.
# * clustering coefficient - measure of the degree to which nodes cluster together.
# * degree centrality - fraction of non-group members connected to group members.
#
# All scores are descriptive statistics of the spatial graph.

###############################################################################
# This dataset contains Leiden cluster groups' annotations in `anndata.AnnData.obs`, which are used for calculation of centrality scores.
#
# First, we need to compute a connectivity matrix from spatial coordinates to calculate the centrality scores. We can use `squidpy.gr.spatial_neighbors` for this purpose. We use the `coord_type="generic"` based on the data and the neighbors are classified with Delaunay triangulation by specifying `delaunay=True`.

sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

###############################################################################
# Compute centrality scores
# -------------------------

###############################################################################
# Centrality scores are calculated with `squidpy.gr.centrality_scores`, with the Leiden groups as clusters.

sq.gr.centrality_scores(adata, cluster_key="leiden")

###############################################################################
# The results were visualized by plotting the average centrality, closeness centrality, and degree centrality using `squidpy.pl.centrality_scores`.

sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(16, 5))

###############################################################################
# ----

###############################################################################
# Compute co-occurrence probability
# ---------------------------------

###############################################################################
# This example shows how to compute the co-occurrence probability.
#
# The co-occurrence score is defined as:
#
# \begin{equation}
# \frac{p(exp|cond)}{p(exp)}
# \end{equation}
# where $p(exp|cond)$ is the conditional probability of observing a
# cluster $exp$ conditioned on the presence of a cluster $cond$, whereas
# $p(exp)$ is the probability of observing $exp$ in the radius size of
# interest. The score is computed across increasing radii size around each
# cell in the tissue.

###############################################################################
# We can compute the co-occurrence score with `squidpy.gr.co_occurrence`.
# Results of co-occurrence probability ratio can be visualized with `squidpy.pl.co_occurrence`. The '3' in the $\frac{p(exp|cond)}{p(exp)}$ represents a Leiden clustered group.

###############################################################################
# We can further visualize tissue organization in spatial coordinates with `squidpy.pl.spatial_scatter`, with an overlay of the expressed genes which were colored in consonance with the Leiden clusters.

adata_subsample = sc.pp.subsample(adata, fraction=0.5, copy=True)

""
sq.gr.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
)
sq.pl.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
    clusters="12",
    figsize=(10, 10),
)
sq.pl.spatial_scatter(
    adata_subsample,
    color="leiden",
    shape=None,
    size=2,
)

###############################################################################
# ----

###############################################################################
# Neighbors enrichment analysis
# -----------------------------

###############################################################################
# This example shows how to run the neighbors enrichment analysis routine.
#
# It calculates an enrichment score based on proximity on the connectivity graph of cell clusters. The number of observed events is compared against $N$ permutations and a *z-score* is computed.

###############################################################################
# This dataset contains cell type annotations in `anndata.Anndata.obs` which are used for calculation of the neighborhood enrichment. We calculate the neighborhood enrichment score with `squidpy.gr.nhood_enrichment`.

sq.gr.nhood_enrichment(adata, cluster_key="leiden")

###############################################################################
# And visualize the results with `squidpy.pl.nhood_enrichment`.

fig, ax = plt.subplots(1, 2, figsize=(13, 7))
sq.pl.nhood_enrichment(
    adata,
    cluster_key="leiden",
    figsize=(8, 8),
    title="Neighborhood enrichment adata",
    ax=ax[0],
)
sq.pl.spatial_scatter(adata_subsample, color="leiden", shape=None, size=2, ax=ax[1])

###############################################################################
# ----

###############################################################################
# Compute Ripley's statistics
# ---------------------------

###############################################################################
# This example shows how to compute the Ripley's L function.
#
# The Ripley's L function is a descriptive statistics function generally used to determine whether points have a random, dispersed or clustered distribution pattern at certain scale. The Ripley's L is a variance-normalized version of the Ripley's K statistic. There are also 2 other Ripley's statistics available (that are closely related): 'G' and 'F'.
#
# Ripley's G monitors the portion of points for which the nearest neighbor is within a given distance threshold, and plots that cumulative percentage against the increasing distance radii.
#
# For increasing separation range, Ripley's F function assembles the percentage of points which can be found in the aforementioned range from an arbitrary point pattern spawned in the expanse of the noticed pattern. 

###############################################################################
# We can compute the Ripley's L function with `squidpy.gr.ripley`.
# Results can be visualized with `squidpy.pl.ripley`. The other Ripley's statistics can be specified using `mode = 'G'` or `mode = 'F'`.

fig, ax = plt.subplots(1, 2, figsize=(15, 7))
mode = "L"

sq.gr.ripley(adata, cluster_key="leiden", mode=mode)
sq.pl.ripley(adata, cluster_key="leiden", mode=mode, ax=ax[0])

sq.pl.spatial_scatter(
    adata_subsample,
    color="leiden",
    groups=["0", "1", "3"],
    shape=None,
    size=2,
    ax=ax[1],
)

###############################################################################
# ----

###############################################################################
# Compute Moran's I score
# -----------------------

###############################################################################
# This example shows how to compute the Moran's I global spatial auto-correlation statistics.
#
# The Moran's I global spatial auto-correlation statistics evaluates whether features (i.e. genes) shows a pattern that is clustered, dispersed or random in the tissue are under consideration.

###############################################################################
# We can compute the Moran's I score with `squidpy.gr.spatial_autocorr` and `mode = 'moran'`. We first need to compute a spatial graph with `squidpy.gr.spatial_neighbors`. We will also subset the number of genes to evaluate.

sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)
sq.gr.spatial_autocorr(
    adata_subsample,
    mode="moran",
    n_perms=100,
    n_jobs=1,
)
adata_subsample.uns["moranI"].head(10)

###############################################################################
# We can visualize some of those genes with `squidpy.pl.spatial_scatter`. We could also pass `mode = 'geary'` to compute a closely related auto-correlation statistic, [Geary's C](https://en.wikipedia.org/wiki/Geary%27s_C). See `squidpy.gr.spatial_autocorr` for more information.

sq.pl.spatial_scatter(
    adata_subsample,
    color=[
        "Slc17a7",
        "Npy2r",
    ],
    shape=None,
    size=2,
    img=False,
)
