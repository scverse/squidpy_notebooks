"""
This is a dummy tutorial
------------------------
"""
from anndata import AnnData
import squidpy as sp

import numpy as np

# %%
# Get the version.
sp.__version__

# %%
# Initialize the :mod:`anndata` object.
r = np.random.RandomState(100)
adata = AnnData(r.rand(200, 100), obs={"cluster": r.randint(0, 3, 200)})  # this is a regular comment
adata

# %%
# This is a new cell. We can reference the docs as :func:`squidpy.gr.moran`.
# Note than any such references or code usage will be automatically linked under that function.
# See example at :ref:`sphx_glr_auto_examples_graph_compute_dummy.py`.
adata.obsm["spatial"] = np.stack([r.randint(0, 500, 200), r.randint(0, 500, 200)], axis=1)
sp.gr.spatial_neighbors(adata, spatial_key="spatial", n_rings=2)
adata
