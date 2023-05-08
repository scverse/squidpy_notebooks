#!/usr/bin/env python
"""
Compute Sepal score
-------------------

This example shows how to compute the Sepal score for spatially variable genes identification.

The Sepal score is a method that simulates a diffusion process to quantify spatial structure in tissue.
See :cite:`andersson2021` for reference.

:::{seealso}

    - See {doc}`sphx_glr_auto_examples_graph_compute_co_occurrence.py` and
      :ref:`sphx_glr_auto_examples_graph_compute_moran.py` for other scores to identify spatially variable genes.
    - See {doc}`sphx_glr_auto_examples_graph_compute_spatial_neighbors.py` for general usage of
      {func}`squidpy.gr.spatial_neighbors`.

:::
"""
import squidpy as sq

adata = sq.datasets.visium_hne_adata()
adata

###############################################################################
# We can compute the Sepal score with {func}`squidpy.gr.sepal`.
# there are 2 important aspects to consider when computing sepal:
#
# - The function only accepts grid-like spatial graphs. Make sure to specify the
#   maximum number of neighbors in your data (6 for an hexagonal grid like Visium)
#   with ``max_neighs = 6``.
# - It is useful to filter out genes that are expressed in very few observations
#   and might be wrongly identified as being spatially variable. If you are performing
#   pre-processing with Scanpy, there is a convenient function that can be used BEFORE
#   normalization {func}`scanpy.pp.calculate_qc_metrics`. It computes several useful
#   summary statistics on both observation and feature axis. We will be using the
#   ``n_cells`` columns in `adata.var` to filter out genes that are expressed in
#   less than 100 observations.
#
# Before computing the Sepal score, we first need to compute a spatial graph with {func}`squidpy.gr.spatial_neighbors`.
# We will also subset the number of genes to evaluate for efficiency purposes.
sq.gr.spatial_neighbors(adata)
genes = adata.var_names[(adata.var.n_cells > 100) & adata.var.highly_variable][0:100]
sq.gr.sepal(adata, max_neighs=6, genes=genes, n_jobs=1)
adata.uns["sepal_score"].head(10)

###############################################################################
# We can visualize some of those genes with {func}`squidpy.pl.spatial_scatter`.
sq.pl.spatial_scatter(adata, color=["Lct", "Ecel1", "Cfap65"])
