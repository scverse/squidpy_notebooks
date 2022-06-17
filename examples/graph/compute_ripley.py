#!/usr/bin/env python
"""
Compute Ripley's statistics
---------------------------

This example shows how to compute the Ripley's L function.

The Ripley's L function is a descriptive statistics generally used
to determine whether points have a random, dispersed or clustered distribution
pattern at certain scale.
The Ripley's L is a variance-normalized version of the Ripley's K statistic.

.. seealso::

    See :ref:`sphx_glr_auto_examples_graph_compute_co_occurrence.py` for
    another score to describe spatial patterns with :func:`squidpy.gr.co_occurrence`.
"""
import squidpy as sq

adata = sq.datasets.slideseqv2()
adata

###############################################################################
# We can compute the Ripley's L function with :func:`squidpy.gr.ripley`.
# Results can be visualized with :func:`squidpy.pl.ripley`.
mode = "L"
sq.gr.ripley(adata, cluster_key="cluster", mode=mode)
sq.pl.ripley(adata, cluster_key="cluster", mode=mode)

###############################################################################
# We can further visualize tissue organization in spatial coordinates
# with :func:`squidpy.pl.spatial_scatter`.
sq.pl.spatial(adata, color="cluster", size=20, shape=None)

###############################################################################
# There are also 2 other Ripley's statistics available (that are closely related):
# ``mode = 'F'`` and ``mode = 'G'``.
