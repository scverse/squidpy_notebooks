#!/usr/bin/env python
"""
Calculate distances to a user-defined anchor point
---------------------------

This example shows how to use :func:`squidpy.tl.exp_dist` to calculate the minimum distances of all observations
to a user-defined anchor point and store the results in :attr:`anndata.AnnData.obsm`.
For plotting expression by distance, see :func:`squidpy.pl.exp_dist`
"""

import squidpy as sq

adata = sq.datasets.mibitof()

###############################################################################
# This data set contains a cell type annotation in :attr:`anndata.AnnData.obs["Cluster"]` and a slide annotation in :attr:`anndata.AnnData.obs["library_id"]`
adata.obs

###############################################################################
# For each slide we now want to calculate the distance of all observations to the closest Epithelial cell.
# In addition we want to include the condition of the donors and the donor id in the resulting design matrix
sq.tl.exp_dist(adata=adata, groups="Epithelial", cluster_key="Cluster", library_key="library_id", covariates=["category","donor"])

###############################################################################
# Since we didn't specify a name, the resulting data frame is stored as "design_matrix" in :attr:`anndata.AnnData.obsm`.
adata.obsm["design_matrix"]

###############################################################################
# NaN values indicate, that the observation belongs to an anchor point. 
# If NaN values are present in the "_raw" distances column as well, then the coordinates for this observation weren't available from the beginning.
