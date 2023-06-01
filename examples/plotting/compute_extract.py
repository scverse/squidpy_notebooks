#!/usr/bin/env python
"""
Plot features in adata.obsm
---------------------------

This example shows how to use :func:`squidpy.pl.extract` to plot features in :attr:`anndata.AnnData.obsm`.

.. seealso::

    See :ref:`sphx_glr_examples_image_compute_summary_features.py` for computing an example of such features.
"""
import squidpy as sq

adata = sq.datasets.slideseqv2()
adata

###############################################################################
# In this dataset, we have saved deconvolution results in :attr:`anndata.AnnData.obsm` and we
# would like to plot them with :func:`squidpy.pl.spatial_scatter`.
adata.obsm["deconvolution_results"].head(10)

###############################################################################
# Squidpy provides an easy wrapper that creates a temporary copy of the
# feature matrix and pass it to :attr:`anndata.AnnData.obs`.
sq.pl.spatial_scatter(
    sq.pl.extract(adata, "deconvolution_results"),
    shape=None,
    color=["Astrocytes", "Mural", "CA1_CA2_CA3_Subiculum"],
    size=4,
)
