#!/usr/bin/env python

# # Visium (spatial statistics and graph)
# This tutorial shows how to apply Squidpy for the analysis of Visium spatial transcriptomics dataset. The dataset used here consist of a Visium slide of a coronal section of the mouse brain. The original dataset is publicly availbale at the 10x genomics `dataset portal https://support.10xgenomics.com/spatial-gene-expression/datasets`_ . Here, we provide a pre-processed dataset, with pre-annoated clusters, in AnnData format. Couple of notes on pre-processing:
# - The pre-processing pipeline is the same as the one showed in the original `Scanpy tutorial https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html`_ .
# - The cluster annotation was performed using several resources, such as the `Allen Brain Atlas http://mouse.brain-map.org/experiment/thumbnails/100048576?image_type=atlas`_ , the `Mouse Brain gene expression atlas(http://mousebrain.org/genesearch.html`_ from the Linnarson lab and this recent `preprint https://www.biorxiv.org/content/10.1101/2020.07.24.219758v1`_ .
#
# ## Import packages & data
# You can run the notebooks in your own conda environemnt create your own with `conda create -f environment.yml`. The file `environment.yml` can be found `here `_ .

import scanpy as sc
import squidpy as sq

import numpy as np

sc.logging.print_header()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3
sq.__version__


# ### Download data
# To get the data, you can simply run the following command, which downloads the relevant data from this [fighare link]() and stores it in `~/.cache`.

adata = sc.read("/Users/giovanni.palla/Datasets/tutorial_data/visium_hne.h5ad")
adata


# ## Analysis and visualization
# First, let's visualize cluster annotation in spatial context with `scanpy.pl.spatial https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.spatial.html`_ .

sc.pl.spatial(adata, color="cluster")


# ### Neighborhood enrichment
# Similar to other spatial data, we can investigate spatial organization of clusters in a quantitative way, by computing a neighborhood enrichment score. You can compute such score with the following function: :func:`squidpy.gr.nhood_enrichment`. In short, it's an enrichment score on spatial proximity of clusters: if spots belonging to two different clusters are often close to each other, then they will have a high score and can be defined as being *enriched*. On the other hand, if they are far apart, and therefore are seldomly neighborhood, the score will be low and they can be defined as *depleted*. This score is based on a permutation-based test, and you can set the number of permutations with the `n_perms` argument (default is 1000).
#
# Since the function works on a connectivity matrix, we need to compute that as well. This can be done with :func:`squidpy.gr.spatial_neighbors`. Please see #ADD LINK for a thorough explanation of how this function works.
#
# Finally, we'll directly visualize the results with :func:`squidpy.pl.nhood_enrichment`.

sq.gr.spatial_neighbors(adata)
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster")


# Given the spatial organization of the mouse brain coronal section, not suprisingly we find that most of the times spots belonging to the same cluster are close to themselves, rather than to other clusters. We do find some enrichment though in the Hippocampus region: *Pyramidal_layer_dentate_gyrus* and *Pyramidal_layer* clusters seems to be often neighbors with the larger *Hippocampus* cluster. This is no news, since that structure has been extensively studied.

# ### Co-occurrence across spatial dimensions
# In addition to the neighbor enrichemnt score, we can visualize cluster co-occurrence in spatial dimensions. This is a simialr analysis of the one presented above, yet it does not operates on the connectivity matrix, yet on the original spatial coordinates. The co-occurrence score is defined as:
# \begin{equation*}
# \frac{p(exp|cond)}{p(exp)}
# \end{equation*}
#
# where $p(exp|cond)$ is the conditional probability of observing a cluster $exp$ conditioned on the presence of a cluster $cond$, whereas $p(exp)$ is the probability of observing $exp$ in the radius size of interest.
#
# We are gonna compute such score with :func:`squidpy.gr.co_occurrence` and set the cluster annotation for the conditional probability with the argument `group`. Then, we visualize the results with :func:`squidpy.pl.co_occurrence`.

sq.gr.co_occurrence(adata, cluster_key="cluster")
sq.pl.co_occurrence(
    adata,
    cluster_key="cluster",
    clusters="Hippocampus",
    figsize=(8, 4),
)


# The result largely recapitulates the previous analysis: the clusters related to the *Pyramidal layers* seem to occu occur at short distances with the larger *Hippocampus* cluster. It should be noted that the distance units are given in pixels of the Visium `source_image`, and corresponds to the spatial coordinates saved in `adata.obsm["spatial"]`

# ### Ligand-receptor interaction analysis
# We are continuing the analysis showing couple of feature-level methods that are very relevant for the analysis of spatial molecular data. For instance, after quantification of cluster co-occurrence, we might be interested in finding molecular instances that could potentially drive such communication. This naturally translates in doing a ligand-receptor interaction analysis. In Squidpy, we provide a fast re-implementation the popular method CellPhoneDB (`paper https://www.nature.com/articles/s41596-020-0292-x`_ - `code https://github.com/Teichlab/cellphonedb`_ ) and extended its databased of annotated ligand-receptor interaction pairs with the popular database `Omnipath https://omnipathdb.org/`_ . You can run the analysis for all clusters pairs, and all genes (in seconds, without leaving this notebook), with :func:`squidpy.gr.ligrec`.
#
# Furthermore, we'll directly visualize the results, filering out low-expressed genes (with the `means_range` argument) and increasing the threshold for the adjusted p-value (with the `alpha` argument). We'll also subset the virualization for only one source group, the *Hippocampus* cluster, and two target groups, *Pyramidal_layer_dentate_gyrus* and *Pyramidal_layer* cluster.

sq.gr.ligrec(
    adata,
    cluster_key="cluster",
)


sq.pl.ligrec(
    adata,
    cluster_key="cluster",
    source_groups="Hippocampus",
    target_groups=["Pyramidal_layer", "Pyramidal_layer_dentate_gyrus"],
    means_range=(3, np.inf),
    alpha=1e-4,
    swap_axes=True,
)


# The dotplot visualization provides an interesting set of candidate ligand-receptor annotation that could be involved in cellular interactions in the hippocampus. A more refined analysis would be for instance to integrate these results with a deconvolution method, to understand what's the proportion of single-cell cell types present in this region of the tissue.

# ### Spatially Variable genes with Moran's I
# Finally, we might be interested in finding genes that show spatial patterns. There are several methods that aimed at address this expliclty, based on point processes or gaussian process regression framework:
# * SPARK `paper https://www.nature.com/articles/s41592-019-0701-7#Abs1`_ - `code https://github.com/xzhoulab/SPARK`_
# * Spatial DE `paper https://www.nature.com/articles/nmeth.4636`_ - `code https://github.com/Teichlab/SpatialDE`_
# * trendsceek `paper https://www.nature.com/articles/nmeth.4634`_ - `code https://github.com/edsgard/trendsceek`_
# * HMRF `paper https://www.nature.com/articles/nbt.4260`_ - `code https://bitbucket.org/qzhudfci/smfishhmrf-py/src/default/`_
#
# Here, we provide a simple approach based on the well-knwon `Moran's I statistics https://en.wikipedia.org/wiki/Moran%27s_I`_ which is in fact used also as a baseline method in the spatially variable gene papers listed above. The function in scanpy is called :func:`squidpy.gr.moran`, and returns both test statistics and adjusted pvalues in `adata.var` slot.

sq.gr.moran(adata, n_jobs=4)


# Thre results are saved in the following `adata.uns` slot. Genes have already been sorted by Moran'I effect size. We can select few genes and visualize their expression levels in the tissue

adata.uns["moranI"].head(10)


sc.pl.spatial(adata, color=["Nrgn", "Camk2n1", "Mobp", "cluster"])


# Interestingly, some of these genes seems to be related to the pyramidal layers, the cortex and the Fiber tract.
