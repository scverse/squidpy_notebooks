# +
"""
This is a test permutations example for four_i example
-----------------------

This is how we can link an example saved under **`examples/graph/compute_non_existing.py`**:
**:ref:`sphx_glr_auto_examples_graph_compute_non_existing.py`**.

All examples should be prefixed with either **compute_** or **plot_**.
"""

import squidpy as sq

# this is for autocompletion and will be removed during python conversion
# # %config Completer.use_jedi = False


adata = sq.datasets.four_i()
adata

# +
# This is a new cell. We can reference the docs as :func:`squidpy.gr.moran`.
# Note than any such references or code usage will be automatically linked under that function.
# See tutorial at :ref:`sphx_glr_auto_tutorials_tutorial_dummy.py`.
# -

sq.gr.spatial_neighbors(adata)

# https://stackoverflow.com/questions/42092218/how-to-add-a-label-to-seaborn-heatmap-color-bar
sq.pl.nhood_enrichment(adata, cluster_key='leiden', figsize=[4, 4], cmap='Blues',
                       cbar_kws={'label': 'Z-score'}, vmin=5, vmax=20)

