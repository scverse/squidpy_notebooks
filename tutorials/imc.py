#!/usr/bin/env python
# coding: utf-8

# # Imaging Mass Cytometry
# This tutorial shows how to apply Squidpy to Imaging Mass Cytometry data. The data used here comes from a recent paper from [Jackson et al.](https://www.nature.com/articles/s41586-019-1876-x). We provide a pre-processed subset of the data, in AnnData format, that can be downloaded [here]().
# 
# ## Import packages & data
# You can run the notebooks in your own conda environemnt create your own with `conda create -f environment.yml`. The file `environment.yml` can be found [here]().

import scanpy as sc
import squidpy as sq

sc.logging.print_header()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3
sq.__version__


# ### Download data
# To get the data, you can simply run the following command, which downloads the relevant data from this [fighare link]() and stores it in `~/.cache`.

adata = sc.read("/Users/giovanni.palla/Datasets/tutorial_data/imc.h5ad")
adata


# ## Analysis and visualization
# First, let's visualize cell type annotation in spatial context. For this, we can use the convenient function `scanpy.pl.spatial`.

sc.pl.spatial(adata, color="cell type", spot_size=10)


# A first glimpse to the data already provides us with a lot of insights. For instance, we can appreciate how the majority of the tissue seems to consist of *apoptotic tumor cells*. There also seems to be other cell types scattered across the tissue of, annotated as *T cells*, *Macrophages* and different types of *Stromal cells*. We can also appreciate how a subset of tumor cell, *basal CK tumor cell* seems to regionalize in clusters in the lower part of the tissue.
# 
# A key goal of Squidpy is to help analyst translates this intuitions on spatial patterns in quantitative scores that provides insight on tissue organization and nrighborhood structure. 

# ### Co-occurrence across spatial dimensions
# In addition to the neighbor enrichemnt score, we can visualize cluster co-occurrence in spatial dimensions. This is a simialr analysis of the one presented above, yet it does not operates on the connectivity matrix, yet on the original spatial coordinates. The co-occurrence score is defined as:
# \begin{equation*}
# \frac{p(exp|cond)}{p(exp)}
# \end{equation*}
# 
# where $p(exp|cond)$ is the conditional probability of observing a cluster $exp$ conditioned on the presence of a cluster $cond$, whereas $p(exp)$ is the probability of observing $exp$ in the radius size of interest. The score is computed across increasing radii size around each cell in the tissue. 
# 
# We are gonna compute such score with :func:`squidpy.gr.co_occurrence` and set the cluster annotation for the conditional probability with the argument `group`. Then, we visualize the results with :func:`squidpy.pl.co_occurrence`. 
# We visualize the result for two conditional groups, namely *basal CK tumor cell* and *T cells*.

sq.gr.co_occurrence(adata, cluster_key="cell type")
sq.pl.co_occurrence(
    adata,
    cluster_key="cell type",
    group=["basal CK tumor cell", "T cells"],
    figsize=(15, 4),
)


# We can observer that *T cells* seems to co-occurr with *endothelial* and *vimentin hi stromal cells*, whereas *basal CK tumor cell* seem to largerly cluster together, except for the presence of a type of stromal cells (*small elongated stromal cell*) at close distance.

# ### Neighborhood enrichment
# A similar analysis that can inform on the neighbor structure of the tissue is the *neighborhood nerichment test*. You can compute such score with the following function: :func:`squidpy.gr.nhood_enrichment`. In short, it's an enrichment score on spatial proximity of clusters: if spots belonging to two different clusters are often close to each other, then they will have a high score and can be defined as being *enriched*. On the other hand, if they are far apart, and therefore are seldomly neighborhood, the score will be low and they can be defined as *depleted*. This score is based on a permutation-based test, and you can set the number of permutations with the `n_perms` argument (default is 1000).
# 
# Since the function works on a connectivity matrix, we need to compute that as well. This can be done with :func:`squidpy.gr.spatial_neighbors`. Please see #ADD LINK for a thorough explanation of how this function works.
# 
# Finally, we'll directly visualize the results with :func:`squidpy.pl.nhood_enrichment`. We will be setting `dendrogram=True` to re-order cell-types rows based on the enrichment score.

sq.gr.spatial_neighbors(adata)
sq.gr.nhood_enrichment(adata, cluster_key="cell type")
sq.pl.nhood_enrichment(adata, cluster_key="cell type")


# Interestingly, *T cells* shows an enrichment with *stromal*  and *endothelial cell*s, as well as *macrophages*. Another interesting result is that *apoptotic tumor cells*, being uniformly spread across the tissue area, show a neighbor depletion against any other cluster (but a strong enrichment for itself). This is a correct interpretation from a permutation based approach, because the cluster annotation, being uniformly spread across the tissue, and in high number, it's more likley to be enriched with cell types from the same class, rather than differen one.
# 
# This is in constrast with a parallel analysis that can be performed with Squidpy, with the function :func:`squidpy.gr.interaction_matrix`. This method builds an interaction matrix based on the number of edges that each cluster share with all the others. By visualizing the results with  :func:`squidpy.pl.interaction_matrix`, we can appreciate that *apoptotic tumor cells* share indeed many edges with other annotations, being uniformly spread across the tissue area. This is a descriptive statistics, but nonetheless useful to understand organization of cell types in tissue.

sq.gr.interaction_matrix(adata, cluster_key="cell type")
sq.pl.interaction_matrix(adata, cluster_key="cell type")


# ### Centrality scores
# Finally, similar to the previous analysis, we can investigate properties of the spatial graph by computing different network centralities network centralities: 
# - degree_centrality
# - clustering_coefficient
# - betweenness_centrality
# - closeness_centrality
# 
# Squidpy provides a convenient function for all of them: :func:`squidpy.gr.centrality_scores` and :func:`squidpy.pl.centrality_scores` for visualization.

sq.gr.centrality_scores(adata, cluster_key="cell type",)
sq.pl.centrality_scores(adata, cluster_key="cell type", figsize=(20, 5), s=500)

