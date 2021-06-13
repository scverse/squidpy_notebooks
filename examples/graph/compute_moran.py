#!/usr/bin/env python
"""
Compute Moran's I score
-----------------------

This example shows how to compute the Moran's I global spatial autocorrelation statistics.

The Moran's I global spatial autocorrelation statistics evaluates whether
features (i.e. genes) shows a pattern that is clustered, dispersed or random
in the tissue are under consideration.

.. seealso::

    - See :ref:`sphx_glr_auto_examples_graph_compute_co_occurrence.py` and
      :ref:`sphx_glr_auto_examples_graph_compute_ripley_k.py` for other scores to describe spatial patterns.
    - See :ref:`sphx_glr_auto_examples_graph_compute_spatial_neighbors.py` for general usage of
      :func:`squidpy.gr.spatial_neighbors`.
"""
import scanpy as sc
import squidpy as sq

adata = sq.datasets.visium_hne_adata()
adata

###############################################################################
# We can compute the Moran's I score with :func:`squidpy.gr.moran`.
# We first need to compute a spatial graph with :func:`squidpy.gr.moran`.
# We will also subset the number of genes to evaluate.

genes = adata[:, adata.var.highly_variable].var_names.values[:100]
sq.gr.spatial_neighbors(adata)
sq.gr.moran(
    adata,
    genes=genes,
    n_perms=100,
    n_jobs=1,
)
adata.uns["moranI"].head(10)

###############################################################################
# We can visualize some of those genes with :func:`scanpy.pl.spatial`.
sc.pl.spatial(adata, color=["Resp18", "Tuba4a"])
