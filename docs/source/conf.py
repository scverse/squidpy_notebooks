# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import os
import sys
from collections import ChainMap
from datetime import datetime
from pathlib import Path

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import squidpy  # TODO
from sphinx_gallery.sorting import ExplicitOrder, _SortKey

sys.path.insert(0, os.path.abspath("_ext"))
needs_sphinx = "3.0"

# -- Project information -----------------------------------------------------

project = "Squidpy"
author = "TODO"  # squidpy.__author__
copyright = f"{datetime.now():%Y}, {author}."
release = "master"
version = "TODO"  # f"master ({squidpy.__version__})"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_gallery.gen_gallery",
    "sphinx_last_updated_by_git",
    "sphinx_copybutton",
    "edit_on_github",
]


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
source_suffix = [".rst", ".ipynb"]
master_doc = "index"
pygments_style = "sphinx"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "**.ipynb",
    "**.md5",
    "**.py",
    "**.ipynb_checkpoints",
]  # ignore anything that isn't .rst

# -- sphinx gallery


def reset_matplotlib(_gallery_conf, _fname):
    import matplotlib as mpl

    mpl.use("agg")

    import matplotlib.pyplot as plt

    plt.rcdefaults()
    mpl.rcParams["savefig.bbox"] = "tight"
    mpl.rcParams["savefig.transparent"] = True


example_dir = Path(__file__).parent.parent.parent / "examples"
tutorial_dir = Path(__file__).parent.parent.parent / "tutorials"
rel_example_dir = Path("..") / ".." / "examples"
rel_tutorial_dir = Path("..") / ".." / "tutorials"


class ExplicitSubsectionOrder(_SortKey):

    _order = ChainMap(
        {
            example_dir / "graph" / "compute_dummy.py": 0,
        },
        {
            tutorial_dir / "tutorial_dummy.py": 0,
        },
    )

    def __call__(self, filename: str) -> int:
        src_file = os.path.normpath(os.path.join(self.src_dir, filename))
        return self._order[Path(src_file)]

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}: {repr(dict(self._order))}>"


sphinx_gallery_conf = {
    "image_scrapers": "matplotlib",  # TODO: napari scraper
    "reset_modules": (
        "seaborn",
        reset_matplotlib,
        # TODO: reset scanpy/napari
    ),
    "filename_pattern": f"{os.path.sep}(plot_|compute_|tutorial_)",
    "examples_dirs": [example_dir, tutorial_dir],
    "gallery_dirs": ["auto_examples", "auto_tutorials"],
    "abort_on_example_error": True,
    "show_memory": True,
    "within_subsection_order": ExplicitSubsectionOrder,
    "subsection_order": ExplicitOrder(
        [
            rel_example_dir / "graph",
            rel_example_dir / "image",
            rel_example_dir / "plot",
        ]
    ),
    "reference_url": {
        "sphinx_gallery": None,
    },
    "line_numbers": False,
    "compress_images": (
        "images",
        "thumbnails",
        "-o3",
    ),  # TODO: CI needs to install optipng
    "inspect_global_variables": False,
    "backreferences_dir": "gen_modules/backreferences",
    "doc_module": "squidpy",
    "download_all_examples": False,
    "pypandoc": True,  # convert rST to md when downloading notebooks
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
todo_include_todos = False

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {"navigation_depth": 4, "logo_only": True}
html_show_sphinx = False

github_repo = "squidpy"
github_repo_nb = "squidpy_notebooks"
