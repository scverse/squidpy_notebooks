#!/usr/bin/env python
"""
Plot features in `adata.obsm`
-----------------------------

This example shows how to use :func:`squidpy.pl.extract` to plot features in `adata.obsm`.

This function is useful when a feature matrix is saved separately in `adata.obsm` and its features
are therefore not accessible via standard scanpy plotting.

.. seealso::

    - See :ref:`sphx_glr_auto_examples_image_compute_summary_features.py`
      for computing an example of such features.
"""
import scanpy as sc
import squidpy as sq

adata = sq.datasets.slideseqv2()
adata

###############################################################################
# In this dataset, we have saved deconvolution results in `adata.obsm` and we
# would like to plot them with :func:`scanpy.pl.plotting`.

adata.obsm["deconvolution_results"].head(10)

###############################################################################
# Squidpy provides an easy wrapper that creates a temporary copy of the
# feature matrix and pass it to `adata.obs` and makes it therefore accessible
# for Scanpy plotting.

sc.pl.spatial(
    sq.pl.extract(adata, "deconvolution_results"),
    color=["Astrocytes", "Mural", "CA1_CA2_CA3_Subiculum"],
    spot_size=30,
)
