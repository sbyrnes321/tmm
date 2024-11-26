# Configuration file for the Sphinx documentation builder.

project = 'tmm'
copyright = '2024, Steven Byrnes'
author = 'Steven Byrnes'

# The full version, including alpha/beta/rc tags
release = '0.1.9'

# Add any Sphinx extension module names here, as strings
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The theme to use for HTML and HTML Help pages.
html_theme = 'sphinx_rtd_theme'

# If true, show URL addresses after external links.
html_show_sourcelink = False
