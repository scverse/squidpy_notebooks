#!/usr/bin/env python
"""
Building spatial neighbors graph
--------------------------------

This example shows how to compute a spatial neighbors graph.

Spatial graph is a graph of spatial neighbors with observations as nodes
and neighbor-hood relations between observations as edges.
We use spatial coordinates of spots/cells to identify neighbors among them.
Different approach of defining a neighborhood relation among observations are used
for different types of spatial datasets.
"""

import scanpy as sc
import squidpy as sq

import numpy as np

###############################################################################
# First, we show how to compute the spatial neighbors graph for a Visium dataset.
adata = sq.datasets.visium_fluo_adata()
adata

###############################################################################
# We use :func:`squidpy.gr.spatial_neighbors` for this.
# The function expects ``coord_type = 'visium'`` by default.
# We set this parameter here explicitly for clarity.
# ``n_rings`` should be used only for Visium datasets.
# It specifies for each spot how many hexagonal rings of spots around
# will be considered neighbors.
sq.gr.spatial_neighbors(adata, n_rings=2, coord_type="grid", n_neighs=6)

###############################################################################
# The function builds a spatial graph and saves its adjacency matrix
# to ``adata.obsp['spatial_connectivities']`` and weighted adjacency matrix to
# ``adata.obsp['spatial_distances']`` by default.
# Note that it can also build a a graph from a square grid, just set ``n_neighs = 4``.
adata.obsp["spatial_connectivities"]

###############################################################################
# The weights of the weighted adjacency matrix are ordinal numbers of hexagonal rings
# in the case of ``coord_type = 'visium'``.
adata.obsp["spatial_distances"]

###############################################################################
# We can visualize the neighbors of a point to better visualize what `n_rings` mean:
_, idx = adata.obsp["spatial_connectivities"][420, :].nonzero()
idx = np.append(idx, 420)
sc.pl.spatial(
    adata[idx, :],
    neighbors_key="spatial_neighbors",
    edges=True,
    edges_width=1,
    img_key=None,
)

###############################################################################
# Next, we show how to compute the spatial neighbors graph for a non-grid dataset.
adata = sq.datasets.imc()
adata

###############################################################################
# We use the same function for this with ``coord_type = 'generic'``.
# ``n_neighs`` and ``radius`` can be used for non-Visium datasets.
# ``n_neighs`` specifies a fixed number of the closest spots for each spot as neighbors.
# Alternatively, ``delaunay = True`` can be used, for a Delaunay triangulation graph.
sq.gr.spatial_neighbors(adata, n_neighs=10, coord_type="generic")
_, idx = adata.obsp["spatial_connectivities"][420, :].nonzero()
idx = np.append(idx, 420)
sc.pl.spatial(
    adata[idx, :],
    color="cell type",
    neighbors_key="spatial_neighbors",
    spot_size=1,
    edges=True,
    edges_width=1,
    img_key=None,
)

###############################################################################
# We use the same function for this with ``coord_type = 'generic'`` and ``delaunay = True``.
# You can appreciate that the neighbor graph is slightly different than before.
sq.gr.spatial_neighbors(adata, delaunay=True, coord_type="generic")
_, idx = adata.obsp["spatial_connectivities"][420, :].nonzero()
idx = np.append(idx, 420)
sc.pl.spatial(
    adata[idx, :],
    color="cell type",
    neighbors_key="spatial_neighbors",
    spot_size=1,
    edges=True,
    edges_width=1,
    img_key=None,
)

###############################################################################
# In order to get all spots within a specified radius (in units of the spatial coordinates)
# from each spot as neighbors, the parameter ``radius`` should be used.
sq.gr.spatial_neighbors(adata, radius=0.3, coord_type="generic")

adata.obsp["spatial_connectivities"]
adata.obsp["spatial_distances"]
