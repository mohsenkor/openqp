# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'OpenQP'
copyright = '2025, OpenQP Team'
author = 'OpenQP Team'
release = 'v1.0-dev'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
import os
import sys
#sys.path.insert(0, os.path.abspath('@CMAKE_SOURCE_DIR@/pyoqp'))  # Adjust path to your Python source
sys.path.insert(0, os.path.abspath('../../pyoqp/oqp'))

extensions = [
    'sphinx.ext.autodoc',    # For auto-generating docs from docstrings
    'sphinx.ext.napoleon',   # For Google/NumPy-style docstrings
]

html_theme = 'sphinx_rtd_theme'
