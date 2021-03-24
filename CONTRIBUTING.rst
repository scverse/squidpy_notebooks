How to write examples and tutorials for squidpy
===============================================

Examples and tutorials are created using sphinx-gallery.
They are automatically run by CI every time a PR is merged in the squidpy repository.

We distinguish between three types of examples and tutorials:

- *Examples* are short explanations of one function (and optionally its related plotting function). See `here <https://squidpy.readthedocs.io/en/latest/auto_examples/graph/compute_interaction_matrix.html>`__ for an example.
- *Tutorials* are longer vignettes, e.g., showing entire workflows. See `here <https://squidpy.readthedocs.io/en/latest/auto_tutorials/tutorial_imc.html>`__ for an example.
- *External tutorials* are tutorials that use external python packages that should be excluded from CI. They should be placed in ``docs/source/external_tutorials/`` and prefixed with ``tutorial_``


Set up environment
------------------

- git clone most recent versions of squidpy and squidpy_notebooks.
- Install latest version of squidpy with ``pip install -e'.[dev,test]'``
- run ``pre-commit install`` in both repos

Datasets for examples/tutorials
-------------------------------

For showcasing functions, please use one of the datasets shipped with squidpy.
Remember that only the visium datasets contain tissue images.

- ``sq.datasets.imc()`` - good for graph
- ``sq.datasets.seqfish()`` - good for graph
- ``sq.datasets.four_i()`` - good for graph
- ``sq.datasets.four_i()`` - good for graph
- ``sq.datasets.visium_fluo_adata()`` - good for graph and image
- ``sq.datasets.visium_hne_adata()`` - good for graph and image
- ``sq.datasets.visium_fluo_adata_crop()`` - good for graph and image (cropped)
- ``sq.datasets.visium_hne_adata_crop()`` - good for graph and image (cropped)
- ``sq.datasets.visium_fluo_image_crop()`` - good for image container (cropped)
- ``sq.datasets.visium_hne_image_crop()`` - good for image container (cropped)
- ``sq.datasets.visium_hne_image()`` - good for image container (cropped)

Main examples and tutorials
---------------------------

Examples and Tutorials are represented as an executable **python file**.
The general structure is described `here <https://sphinx-gallery.github.io/stable/syntax.html>`_ .
You can work on a jupyter notebook to develop the example, but the file needs to be pushed as a ``.py`` file.
You can conveniently go back and forth with ``jupytext`` (install it with pip):

.. code-block::

   jupytext --to py:sphinx mynotebook.ipynb  # create .py file
   jupytext --to notebook mynotebook.py  # create .ipynb file

The pre-commit will flag potential problems for you.
One of the most common problems is that lines of text go over 120 characters, so make sure to check that.
Remember that the text cells in examples will be rendered with rst, so checkout this
`cheatsheet <https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst>`_
You can locally generate the docs to check that everything looks good: ``tox -e docs``.

Make sure to follow the following checklist before merging a new example / tutorial:

- if using math expressions, ensure they render properly (e.g. using the ``:math:`` directive for rst)
- make sure we're referring to the package always the same, e.g. *Squidpy*
- use the ``.. seealso::`` directive to highlight the prominence of other examples in the introductory text
- ensure examples/tutorials are properly linked (sphinx will throw warnings if not).
  Link to examples using the following syntax ``:ref:`sphx_glr_auto_tutorials_tutorial_seqfish.py```.
- ensure that in ``.py`` files, first line after the title is of the following format::

    """
    Super tutorial title
    --------------------
    This {example,tutorial} shows how to ...

    Keep 1 line above free; here comes the general description.
    ... seealso::
    <here recommend other examples>
    """

  The first line after title should be short, since this is the hover info displayed when hovering over the tutorial
- ensure that citations are in ``docs/source/references.bib`` and are used within the examples/tutorials.
  Cite in rst using the ``:cite:`` directive.
- ensure when referencing functions/classes/packages/etc., we use RST roles, such as:
  ``:func:squidpy.im.process_img``, ``:class:squidpy.im.ImageContainer``, etc.
- ensure that example/tutorial titles are capitalized, but do not follow Camel Case style
  (i.e. Process image is good, Process Image is bad)
- for example values, use ``foo = 'bar'`` instead of ``foo='bar'`` or ``foo = "bar"``
  (to be more consistent with main docs)
- ensure that .py examples/tutorials are executable (``chmod +x``) and have ``#!/usr/bin/env python`` shebang at the top


External tutorials
------------------

External tutorials are not run by CI, and therefore provided as **ipynb notebooks**.
Therefore, they will be rendered as they appear locally.

You can link to other `notebooks <https://nbsphinx.readthedocs.io/en/0.8.1/markdown-cells.html#Links-to-Other-Notebooks>`__,
or to other `.rst files <https://nbsphinx.readthedocs.io/en/0.8.1/markdown-cells.html#Links-to-*.rst-Files-(and-Other-Sphinx-Source-Files)>`__
or `functions <https://nbsphinx.readthedocs.io/en/0.8.1/markdown-cells.html#Links-to-Domain-Objects>`__ within them as well (but it's more brittle and inconvenient).

Make sure the tutorial has a thumbnail, see `here <https://nbsphinx.readthedocs.io/en/dask-theme/gallery/cell-metadata.html>`__.
You can also define the hover info.
Put this in the metadata of the cell that produces the image that should be used as thubnail image

.. code-block::

  {
      "nbsphinx-thumbnail": {
          "tooltip": "This tooltip message was defined in cell metadata"
      }
  }
