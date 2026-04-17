# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from datetime import datetime
import re
import sys
import os

# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

sys.path.insert(0, os.path.abspath("."))
sys.path.insert(0, os.path.abspath("../"))
sys.path.insert(0, os.path.abspath("./source/"))


# -- Project information -----------------------------------------------------

from flat import __version__

project = "FlatCAT"
author = "Jure Cerar"
copyright = f"2023-{datetime.now().year}, " + author
packageversion = __version__
release = packageversion

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx_rtd_theme",
]

templates_path = ["templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autodoc_mock_imports = ["numpy", "pymol"]

napoleon_google_docstring = True
napoleon_numpy_docstring = True  # optional, if you're not using NumPy style

napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
html_static_path = ["static"]

html_theme = "sphinx_rtd_theme"
# pygments_style = "sphinx"

def process_signature(app, what, name, obj, options, signature, return_annotation):
    if signature:
        signature = re.sub(
            r"_self=<module 'pymol\.cmd' from '[^']+'>", "_self=cmd", signature)
    return signature, return_annotation


def setup(app):
    app.connect("autodoc-process-signature", process_signature)
