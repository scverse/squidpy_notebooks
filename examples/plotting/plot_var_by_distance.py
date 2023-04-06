#!/usr/bin/env python
"""
Calculate distances to a user-defined anchor point
---------------------------

This example shows how to use :func:`squidpy.pl.var_by_distance` to plot gene expression by distance
to a user-defined anchor point, which is computed by :func:`squidpy.tl.var_by_distance`.
"""

import squidpy as sq

###############################################################################
# Load the data set
adata = sq.datasets.mibitof()

# compute the distances to the closest Epithelial cell for each observation in the data set and include some covariates.
sq.tl.var_by_distance(
    adata=adata, groups="Epithelial", cluster_key="Cluster", library_key="library_id", covariates=["category", "donor"]
)

###############################################################################
# Plot the expression of CD98 by distance to the closest Epithelial cell for each donor.

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
