#!/usr/bin/env python
"""
Analyze Slide-seqV2 data
========================

This tutorial shows how to apply Squidpy for the analysis of Slide-seqV2 data.

The data used here was obtained from :cite:`Stickels2020-rf`.
We provide a pre-processed subset of the data, in :class:`anndata.AnnData` format.
We would like to thank @tudaga for providing cell-type level annotation.
For details on how it was pre-processed, please refer to the original paper.

.. seealso::

    See :ref:`sphx_glr_auto_tutorials_tutorial_visium_hne.py` and
    :ref:`sphx_glr_auto_tutorials_tutorial_seqfish.py` for additional analysis examples.

Import packages & data
----------------------
To run the notebook locally, create a conda environment as *conda env create -f environment.yml* using this
`environment.yml <https://github.com/scverse/squidpy_notebooks/blob/main/environment.yml>`_.
"""

import scanpy as sc
import squidpy as sq

sc.logging.print_header()
print(f"squidpy=={sq.__version__}")

# load the pre-processed dataset
adata = sq.datasets.slideseqv2()
adata

###############################################################################
# First, let's visualize cluster annotation in spatial context
# with :func:`scanpy.pl.spatial`.
sc.pl.spatial(adata, color="cluster", spot_size=30)

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
# Please see :ref:`sphx_glr_auto_examples_graph_compute_spatial_neighbors.py` and
# :ref:`sphx_glr_auto_examples_graph_compute_nhood_enrichment.py` for more details
# of how these functions works.
#
# Finally, we'll directly visualize the results with :func:`squidpy.pl.nhood_enrichment`.
# We'll add a dendrogram to the heatmap computed with linkage method *ward*.

sq.gr.spatial_neighbors(adata, coord_type="generic")
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster", method="single", cmap="inferno", vmin=-50, vmax=100)

###############################################################################
# Interestingly, there seems to be an enrichment between the *Endothelial_Tip*,
# the *Ependymal* cells. Another putative enrichment is between the *Oligodendrocytes*
# and *Polydendrocytes* cells. We can visualize the spatial organization of such clusters.
# For this, we'll use :func:`scanpy.pl.spatial` again.

sc.pl.spatial(
    adata,
    color="cluster",
    groups=["Endothelial_Tip", "Ependymal", "Oligodendrocytes", "Polydendrocytes"],
    spot_size=30,
)

###############################################################################
# Ripley's statistics
# -------------------
# In addition to the neighbor enrichment score, we can further investigate spatial
# organization of cell types in tissue by means of the Ripley's statistics.
# Ripley's statistics allow analyst to evaluate whether a discrete annotation (e.g. cell-type)
# appears to be clustered, dispersed or randomly distributed on the area of interest.
# In Squidpy, we implement three closely related Ripley's statistics, that can be
# easily computed with :func:`squidpy.gr.ripley`. Here, we'll showcase the Ripley's L statistic,
# which is a variance-stabilized version of the Ripley's K statistics.
# We'll visualize the results with :func:`squidpy.pl.ripley`.
# Check :ref:`sphx_glr_auto_examples_graph_compute_ripley.py` for more details.
mode = "L"
sq.gr.ripley(adata, cluster_key="cluster", mode=mode, max_dist=500)
sq.pl.ripley(adata, cluster_key="cluster", mode=mode)

###############################################################################
# The plot highlight how some cell-types have a more clustered pattern,
# like *Astrocytes* and *CA11_CA2_CA3_Subiculum* cells, whereas other have a more
# dispersed pattern, like *Mural* cells. To confirm such interpretation, we can
# selectively visualize again their spatial organization.
sc.pl.spatial(
    adata,
    color="cluster",
    groups=["Mural", "CA1_CA2_CA3_Subiculum"],
    spot_size=30,
)

###############################################################################
# Ligand-receptor interaction analysis
# ------------------------------------
# The analysis showed above has provided us with quantitative information on
# cellular organization and communication at the tissue level.
# We might be interested in getting a list of potential candidates that might be driving
# such cellular communication.
# This naturally translates in doing a ligand-receptor interaction analysis.
# In Squidpy, we provide a fast re-implementation the popular method CellPhoneDB :cite:`cellphonedb`
# (`code <https://github.com/Teichlab/cellphonedb>`_ )
# and extended its database of annotated ligand-receptor interaction pairs with
# the popular database *Omnipath* :cite:`omnipath`.
# You can run the analysis for all clusters pairs, and all genes (in seconds,
# without leaving this notebook), with :func:`squidpy.gr.ligrec`.
#
# Let's perform the analysis and visualize the result for three clusters of
# interest: *Polydendrocytes* and *Oligodendrocytes*.
# For the visualization, we will filter out annotations
# with low-expressed genes (with the ``means_range`` argument)
# and decreasing the threshold
# for the adjusted p-value (with the ``alpha`` argument)
# Check :ref:`sphx_glr_auto_examples_graph_compute_ligrec.py` for more details.
sq.gr.ligrec(
    adata,
    n_perms=100,
    cluster_key="cluster",
    clusters=["Polydendrocytes", "Oligodendrocytes"],
)
sq.pl.ligrec(
    adata,
    cluster_key="cluster",
    source_groups="Oligodendrocytes",
    target_groups=["Polydendrocytes"],
    pvalue_threshold=0.05,
    swap_axes=True,
)

###############################################################################
# The dotplot visualization provides an interesting set of candidate interactions
# that could be involved in the tissue organization of the cell types of interest.
# It should be noted that this method is a pure re-implementation of the original
# permutation-based test, and therefore retains all its caveats
# and should be interpreted accordingly.

###############################################################################
# Spatially variable genes with spatial autocorrelation statistics
# ----------------------------------------------------------------
# Lastly, with Squidpy we can investigate spatial variability of gene expression.
# :func:`squidpy.gr.spatial_autocorr` conveniently wraps two
# spatial autocorrelation statistics: *Moran's I* and *Geary's C**.
# They provide a score on the degree of spatial variability of gene expression.
# The statistic as well as the p-value are computed for each gene, and FDR correction
# is performed. For the purpose of this tutorial, let's compute the *Moran's I* score.
# See :ref:`sphx_glr_auto_examples_graph_compute_moran.py` for more details.
sq.gr.spatial_autocorr(adata, mode="moran")
adata.uns["moranI"].head(10)

###############################################################################
# The results are stored in `adata.uns["moranI"]` and we can visualize selected genes
# with :func:`scanpy.pl.spatial`.
sc.pl.spatial(
    adata,
    color=["Ttr", "Plp1", "Mbp", "Hpca", "Enpp2"],
    spot_size=30,
)
