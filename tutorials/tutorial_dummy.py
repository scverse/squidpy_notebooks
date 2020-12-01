"""
This is a dummy tutorial
------------------------
"""
import numpy as np
from anndata import AnnData

# import squidpy as sp  # TODO: import me once it's public

# %%
# Initialize the :mod:`anndata` object.
r = np.random.RandomState(100)
adata = AnnData(r.rand(200, 100), obs={"cluster": r.randint(0, 3, 200)})  # this is a regular comment
adata

# %%
# This is a new cell. We can reference the docs as :func:`squidpy.graph.moran`.
# Note than any such references or code usage will be automatically linked under that function.
# See example at :ref:`sphx_glr_auto_examples_graph_compute_dummy.py`.
adata.obsm["spatial"] = np.stack([r.randint(0, 500, 200), r.randint(0, 500, 200)], axis=1)
# sp.graph.spatial_connectivity(adata, obsm="spatial", n_rings=2)  # TODO
adata
