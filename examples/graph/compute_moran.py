#!/usr/bin/env python
"""
Compute Moran's I score
-----------------------

This example shows how to compute the Moran's I global spatial auto-correlation statistics.

The Moran's I global spatial auto-correlation statistics evaluates whether
features (i.e. genes) shows a pattern that is clustered, dispersed or random
in the tissue are under consideration.

.. seealso::

    - See :ref:`sphx_glr_examples_graph_compute_co_occurrence.py` and
      :ref:`sphx_glr_examples_graph_compute_ripley.py` for other scores to describe spatial patterns.
    - See :ref:`sphx_glr_examples_graph_compute_spatial_neighbors.py` for general usage of
      :func:`squidpy.gr.spatial_neighbors`.
"""
import squidpy as sq

adata = sq.datasets.visium_hne_adata()
adata

###############################################################################
# We can compute the Moran's I score with :func:`squidpy.gr.spatial_autocorr` and ``mode = 'moran'``.
# We first need to compute a spatial graph with :func:`squidpy.gr.spatial_neighbors`.
# We will also subset the number of genes to evaluate.
genes = adata[:, adata.var.highly_variable].var_names.values[:100]
sq.gr.spatial_neighbors(adata)
sq.gr.spatial_autocorr(
    adata,
    mode="moran",
    genes=genes,
    n_perms=100,
    n_jobs=1,
)
adata.uns["moranI"].head(10)

###############################################################################
# We can visualize some of those genes with :func:`squidpy.pl.spatial_scatter`.
sq.pl.spatial_scatter(adata, color=["Resp18", "Tuba4a"])

###############################################################################
# We could've also passed ``mode = 'geary'`` to compute a closely related auto-correlation statistic, `Geary's C
# <https://en.wikipedia.org/wiki/Geary%27s_C>`_. See :func:`squidpy.gr.spatial_autocorr` for more information.
