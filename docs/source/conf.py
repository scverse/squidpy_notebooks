# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import os
import sys
from datetime import datetime
from pathlib import Path

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import squidpy
from sphinx.application import Sphinx

sys.path.insert(0, os.path.abspath("_ext"))
needs_sphinx = "3.0"

# -- Project information -----------------------------------------------------

project = "Squidpy"
author = squidpy.__author__
copyright = f"{datetime.now():%Y}, {author}."
release = "master"
version = f"{release} ({squidpy.__version__})"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "nbsphinx",
    "sphinx.ext.intersphinx",
    "sphinx_gallery.gen_gallery",
    "sphinx_last_updated_by_git",
    "sphinxcontrib.bibtex",
    "sphinx_copybutton",
    "edit_on_github",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    # TODO: uncomment once the docs are up
    # "squidpy": ("https://squidpy.readthedocs.io/en/stable/", None),
    # TODO: remove me if not used
    "anndata": ("https://anndata.readthedocs.io/en/stable/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable/", None),
    "napari": ("https://napari.org/docs/dev/", None),
    "skimage": ("https://scikit-image.org/docs/stable/", None),
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
source_suffix = [".rst", ".ipynb"]
master_doc = "index"
pygments_style = "sphinx"

# citation
bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = ["author_year"]
bibtex_default_style = "alpha"

# spelling
spelling_lang = "en_US"
spelling_warning = True
spelling_word_list_filename = "spelling_wordlist.txt"
spelling_add_pypi_package_names = True
spelling_show_suggestions = True
spelling_exclude_patterns = ["references.rst"]
# see: https://pyenchant.github.io/pyenchant/api/enchant.tokenize.html
spelling_filters = ["enchant.tokenize.URLFilter", "enchant.tokenize.EmailFilter"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "auto_*/**.ipynb",
    "auto_*/**.md5",
    "auto_*/**.py",
    "**.ipynb_checkpoints",
]  # ignore anything that isn't .rst or .ipynb

# -- sphinx gallery


def reset_matplotlib(_gallery_conf, _fname):
    import matplotlib as mpl

    mpl.use("agg")

    import matplotlib.pyplot as plt

    plt.rcdefaults()
    mpl.rcParams["savefig.bbox"] = "tight"
    mpl.rcParams["savefig.transparent"] = True
    mpl.rcParams["figure.figsize"] = (12, 8)
    mpl.rcParams["figure.dpi"] = 90
    mpl.rcParams["figure.autolayout"] = True


_root = Path(__file__).parent.parent.parent
sphinx_gallery_conf = {
    "image_scrapers": "matplotlib",
    "reset_modules": (
        "seaborn",
        reset_matplotlib,
    ),
    "filename_pattern": f"{os.path.sep}(compute_|tutorial_)",
    "examples_dirs": [_root / "examples", _root / "tutorials"],
    "gallery_dirs": ["auto_examples", "auto_tutorials"],
    "abort_on_example_error": True,
    "show_memory": True,
    "reference_url": {
        "sphinx_gallery": None,
    },
    "line_numbers": False,
    "compress_images": (
        "images",
        "thumbnails",
        "-o3",
    ),
    "remove_config_comments": True,
    "inspect_global_variables": False,
    "backreferences_dir": "gen_modules/backreferences",
    "doc_module": "squidpy",
    "download_all_examples": False,
    "show_signature": False,
    "pypandoc": {
        "extra_args": [
            "--mathjax",
        ],
        "filters": [str(_root / ".scripts" / "filters" / "strip_interpreted_text.py")],
    },
}
nbsphinx_thumbnails = {
    "auto_**": "_static/img/squidpy_vertical.png",
    "external_tutorials/**": "_static/img/squidpy_vertical.png",
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


def setup(app: Sphinx) -> None:
    app.add_css_file("css/custom.css")
