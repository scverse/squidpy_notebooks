#!/usr/bin/env python
"""
Compute Ripley K function
---------------------------

This example shows how to compute the Ripley's K function.

The Ripley K function is a descriptive statistics that is generally used
to determine whether points have a random, dispersed or clustered distribution
pattern at certain scale.

.. seealso::

    See :ref:`sphx_glr_auto_examples_graph_compute_co_occurrence.py` for
    another score to describe spatial patterns with :func:`squidpy.gr.co_occurrence`.
"""
import scanpy as sc
import squidpy as sq

adata = sq.datasets.imc()
adata

###############################################################################
# We can compute the Ripley K function with :func:`squidpy.gr.ripley_k`.
# Results can be visualized with :func:`squidpy.pl.ripley_k`.
sq.gr.ripley_k(adata, cluster_key="cell type")
sq.pl.ripley_k(adata, cluster_key="cell type")

###############################################################################
# We can further visualize tissue organization in spatial coordinates
# with :func:`scanpy.pl.spatial`.
sc.pl.spatial(adata, color="cell type", spot_size=10)
