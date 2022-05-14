#!/usr/bin/env python
"""
Analyze 4i data
===============

This tutorial shows how to apply Squidpy for the analysis of 4i data.

The data used here was obtained from :cite:`Gut2018-ap`.
We provide a pre-processed subset of the data, in :class:`anndata.AnnData` format.
For details on how it was pre-processed, please refer to the original paper.

.. seealso::

    See :ref:`sphx_glr_auto_tutorials_tutorial_imc.py` for additional analysis examples.

Import packages & data
----------------------
To run the notebook locally, create a conda environment as *conda env create -f environment.yml* using this
`environment.yml <https://github.com/theislab/squidpy_notebooks/blob/main/environment.yml>`_.
"""

import scanpy as sc
import squidpy as sq

sc.logging.print_header()
print(f"squidpy=={sq.__version__}")

# load the pre-processed dataset
adata = sq.datasets.four_i()

###############################################################################
# First, let's visualize cluster annotation in spatial context
# with :func:`scanpy.pl.spatial`.
sc.pl.spatial(adata, color="cluster", spot_size=1)


###############################################################################
# Neighborhood enrichment analysis
# --------------------------------
# Similar to other spatial data, we can investigate spatial organization of clusters
# in a quantitative way, by computing a neighborhood enrichment score.
# You can compute such score with the following function: :func:`squidpy.gr.nhood_enrichment`.
# In short, it's an enrichment score on spatial proximity of clusters:
# if spots belonging to two different clusters are often close to each other,
# then they will have a high score and can be defined as being *enriched*.
# On the other hand, if they are far apart, the score will be low
# and they can be defined as *depleted*.
# This score is based on a permutation-based test, and you can set
# the number of permutations with the `n_perms` argument (default is 1000).
#
# Since the function works on a connectivity matrix, we need to compute that as well.
# This can be done with :func:`squidpy.gr.spatial_neighbors`.
# Please see :ref:`sphx_glr_auto_examples_graph_compute_spatial_neighbors.py` for more details
# of how this function works.
#
# Finally, we'll directly visualize the results with :func:`squidpy.pl.nhood_enrichment`.
# We'll add a dendrogram to the heatmap computed with linkage method *ward*.
sq.gr.spatial_neighbors(adata, coord_type="generic")
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster", method="ward", vmin=-100, vmax=100)


###############################################################################
# A similar analysis can be performed with :func:`squidpy.gr.interaction_matrix`.
# The function computes the number of shared edges in the neighbor graph between clusters.
# Please see :ref:`sphx_glr_auto_examples_graph_compute_interaction_matrix.py` for more details
# of how this function works.
sq.gr.interaction_matrix(adata, cluster_key="cluster")
sq.pl.interaction_matrix(adata, cluster_key="cluster", method="ward", vmax=20000)


###############################################################################
# Additional analyses to gain quantitative understanding of spatial patterning of
# sub-cellular observations are:
# - :ref:`sphx_glr_auto_examples_graph_compute_ripley.py` for Ripley's statistics.
# - :ref:`sphx_glr_auto_examples_graph_compute_co_occurrence.py` for co-occurrence score.
#
# Spatially variable genes with spatial autocorrelation statistics
# ----------------------------------------------------------------
# With Squidpy we can investigate spatial variability of gene expression.
# This is an example of a function that only supports 2D data.
# :func:`squidpy.gr.spatial_autocorr` conveniently wraps two
# spatial autocorrelation statistics: *Moran's I* and *Geary's C*.
# They provide a score on the degree of spatial variability of gene expression.
# The statistic as well as the p-value are computed for each gene, and FDR correction
# is performed. For the purpose of this tutorial, let's compute the *Moran's I* score.
# See :ref:`sphx_glr_auto_examples_graph_compute_moran.py` for more details.

adata.var_names_make_unique()
sq.gr.spatial_autocorr(adata, mode="moran")
adata.uns["moranI"].head(10)

###############################################################################
# The results are stored in `adata.uns['moranI']` and we can visualize selected genes
# with :func:`scanpy.pl.spatial`.

sc.pl.spatial(adata, color="Yap/Taz", spot_size=1)
