"""
Receptor-ligand analysis
------------------------
This example shows how to run the receptor-ligand analysis.

It uses an efficient re-implementation of the [CellPhoneDB20]_ algorithm which can handle large
number of interacting pairs (100k+) and cluster combinations (100+).
"""
import squidpy as sq

adata = sq.datasets.seqfish()
adata

# %%
# To get started, we just need an :class:`anndata.AnnData` object with some clustering information. Below are some
# useful parameters of :func:`squidpy.gr.ligrec`:
#
# - ``n_perms`` - number of permutations for the permutation test.
# - ``interactions`` - list of interaction, by default we fetch all available interactions from [OmniPath16]_.
# - ``{interactions,transmitter,receiver}_params`` - parameters used if downloading the ``interactions``,
#   see :func:`omnipah.interactions.import_intercell_network` for more information.
# - ``threshold`` - percentage of cells required to be expressed in a given cluster.
# - ``corr_method`` - false discovery rate (FDR) correction method to use.
#
# Since we're interested in receptors and ligands in this example, we specify these categories in ``receiver_params``
# and ``transmitter_params``, respectively.
# If desired, we can also restrict the resources to just a select few. For example, in order to only use
# [CellPhoneDB20]_, set ``interactions_params={'resources': 'CellPhoneDB'}``.
#
res = sq.gr.ligrec(
    adata,
    n_perms=1000,
    cluster_key="celltype_mapped_refined",
    copy=True,
    use_raw=False,
    transmitter_params={"categories": "ligand"},
    receiver_params={"categories": "receptor"},
)

# %%
# First, we inspect the calculated means. The resulting object is a :class:`pandas.DataFrame`, with rows corresponding
# to interacting pairs and columns to cluster combinations.
res.means.head()

# %%
# Next, we take a look at the p-values. If ``corr_method != None``, this will contained the corrected p-values.
# The p-values marked as `NaN` correspond to interactions, which did not pass the filtering ``threshold`` specified
# above.
res.pvalues.head()

# %%
# Any interaction metadata downloaded from :mod:`omnipath`, such as the interaction type, can be accessed as:
res.metadata.head()

# %%
# In order to plot the results, we can run :func:`squidpy.pl.ligrec`. Some useful parameters are:
#
# - ``{source,target}_groups`` - only plot specific source/target clusters.
# - ``dendrogram`` - whether to hierarchically cluster the rows, columns or both.
# - ``mean_range`` - plot only interactions whose means are in this range.
# - ``pval_threshold`` - plot only interactions whose p-values are below this threshold.
#
# In the plot below, to highlight significance, we've marked all p-values <= 0.005 with tori.
sq.pl.ligrec(res, source_groups="Erythroid", alpha=0.005)
