"""
Neighbors enrichment analysis
------------------------
This example shows how to run the neighbors enrichment analysis routine in squidpy.
It calculates based on pre-defined clusters the proximity between those in the calculated connectivity graph.
The number of observed events is compared versus permutations and Z-scores are summarized.
"""

import squidpy as sq

adata = sq.datasets.visium_fluo_adata()
adata

# %%
# This dataset contains cell type annotations in adata.obs
# which are used for calculation of the neigborhood enrichment.
# First, we need to compute a connectivity matrix from spatial coordinates.

sq.gr.spatial_neighbors(adata)

# %%
# Then we can calculate the neighborhood enrichment score with :func:`squidpy.gr.nhood_enrichment`

sq.gr.nhood_enrichment(adata, cluster_key="cluster")

# %%
# And visualize the results with :func:`squidpy.pl.nhood_enrichment`

sq.pl.nhood_enrichment(adata, cluster_key="cluster")
