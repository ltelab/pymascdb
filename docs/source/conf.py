# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))




# -- Project information -----------------------------------------------------

project = 'mascdb'
copyright = '2021, GG & JG'
author = 'Gionata Ghiggi, Jacopo Grazioli'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
#extensions = ['sphinx.ext.napoleon','sphinx.ext.autosummary']
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
    "sphinx_gallery.gen_gallery",
]


#, 'sphinx.ext.coverage', 'sphinx.ext.autodoc']
#master_doc = 'index'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# --- Mock imports -------------------------------------------------------------
import mock 
 

MOCK_MODULES = ['numpy','scipy','matplotlib', 'matplotlib.pyplot', 'pandas',
 		'dask','xarray','seaborn','dask.diagnostics','skimage',
 		'skimage.morphology','skimage.morphology','skimage.filters']
for mod_name in MOCK_MODULES:
     sys.modules[mod_name] = mock.MagicMock()



# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
