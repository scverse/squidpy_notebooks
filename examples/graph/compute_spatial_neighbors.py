"""
Building spatial neighbors graph
------------------------
This example shows how to compute a spatial neighbors graph.

Spatial graph is a graph of spatial neighbors with spots as nodes
and neighbor-hood relations between spots as edges.
We use spatial coordinates of spots to identify neighborsamong them.
Different approach of defining a neighborhood relation among spots are used
for different types of spatial datasets.
"""

import squidpy as sq

###############################################################################
# First, we show how to compute the spatial neighbors garph for a visium dataset.

adata = sq.datasets.visium_fluo_adata()
adata

###############################################################################
# We use :func:`squidpy.gr.spatial_neighbors` for this.
# The function expects ``coord_type="visium"`` by default.
# We set this parameter here explicitly for clarity.
# ``n_rings`` should be used only for visium datasets.
# It specifies for each spot how many hexagonal rings of spots around
# it will be considered neighbors.

sq.gr.spatial_neighbors(adata, n_rings=2, coord_type="visium")

###############################################################################
# The function builds a spatial graph and saves its adjacency matrix
# to ``adata.obsp["spatial_connectivities"]`` and weighted adjacency matrix to
# ``adata.obsp["spatial_distances"]`` by default.

adata.obsp["spatial_connectivities"]

###############################################################################
# For ``n_rings=1`` there will be no ``adata.obsp["spatial_distances"]``
# The weights of the weighted adjacency matrix are ordinal numbers of hexagonal rings
# in the case of ``coord_type="visium"``.

adata.obsp["spatial_distances"]

###############################################################################
# Next, we show how to compute the spatial neighbors garph for a non-visium dataset.

adata = sq.datasets.visium_fluo_adata()
adata

###############################################################################
# We use the same function for this with ``coord_type="generic"``.
# ``n_neigh`` and ``radius`` can be used for non-visium datasets.
# ``n_neigh`` specifies a fixed number of the closest spots for each spot as neighbors.

sq.gr.spatial_neighbors(adata, n_neigh=10, coord_type="generic")

adata.obsp["spatial_connectivities"]
adata.obsp["spatial_distances"]

###############################################################################
# In order to get all spots within a specified radius (in units of the spatial coordinates)
# from each spot as neighbors, the parameter ``radius`` should be used.

sq.gr.spatial_neighbors(adata, radius=0.3, coord_type="generic")

adata.obsp["spatial_connectivities"]
adata.obsp["spatial_distances"]
