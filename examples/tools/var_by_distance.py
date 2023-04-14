#!/usr/bin/env python
"""
Calculate distances to a user-defined anchor point
---------------------------

This example shows how to use :func:`squidpy.tl.exp_dist` to calculate the minimum distances of all observations
to a user-defined anchor point, store the results in :attr:`anndata.AnnData.obsm` and plot the expression by distance.
using :func:`squidpy.pl.exp_dist`.
"""

import squidpy as sq

adata = sq.datasets.mibitof()

###############################################################################
# This data set contains a cell type annotation in :attr:`anndata.AnnData.obs["Cluster"]`
# and a slide annotation in :attr:`anndata.AnnData.obs["library_id"]`
adata.obs

###############################################################################
# For each slide we now want to calculate the distance of all observations to the closest Epithelial cell.
# In addition we want to include the condition of the donors and the donor id in the resulting design matrix
# As we don't create a copy, the result will be stored in :attr:`anndata.AnnData.obsm`.
sq.tl.exp_dist(
    adata=adata, groups="Epithelial", cluster_key="Cluster", library_key="library_id", covariates=["category", "donor"]
)

###############################################################################
# Since we didn't specify a name, the resulting data frame is called "design_matrix".
# NaN values indicate, that the observation belongs to an anchor point
# or that the coordinates for this observation weren't available from the beginning on.
adata.obsm["design_matrix"]

###############################################################################
# We now want to visualize the results and plot the expression of CD98 by distance to the closest Epithelial cell.
# In addition we want to differentiate between two expression trends by specifying a covariate.
sq.pl.var_by_distance(
    adata=adata,
    design_matrix_key="design_matrix",
    var="CD98",
    anchor_key="Epithelial",
    covariate="donor",
    line_palette=["blue", "orange"],
    show_scatter=False,
    figsize=(5, 4),
)
