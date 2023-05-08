#!/usr/bin/env python
r"""
Compute co-occurrence probability
---------------------------------

This example shows how to compute the co-occurrence probability.

The co-occurrence score is defined as:

.. math::
    \frac{p(exp|cond)}{p(exp)}

where :math:`p(exp|cond)` is the conditional probability of observing a cluster :math:`exp` conditioned
on the presence
of a cluster :math:`cond`, whereas :math:`p(exp)` is the probability of observing :math:`exp` in the
radius size of interest.
The score is computed across increasing radii size around each cell in the tissue.

:::{seealso}

    See {doc}`sphx_glr_auto_examples_graph_compute_ripley.py` for
    another score to describe spatial patterns with {func}`squidpy.gr.ripley`.

:::
"""
import squidpy as sq

adata = sq.datasets.imc()
adata

###############################################################################
# We can compute the co-occurrence score with {func}`squidpy.gr.co_occurrence`.
# Results can be visualized with {func}`squidpy.pl.co_occurrence`.
sq.gr.co_occurrence(adata, cluster_key="cell type")
sq.pl.co_occurrence(adata, cluster_key="cell type", clusters="basal CK tumor cell")

###############################################################################
# We can further visualize tissue organization in spatial coordinates with {func}`squidpy.pl.spatial_scatter`.
sq.pl.spatial_scatter(adata, color="cell type", size=10, shape=None)
