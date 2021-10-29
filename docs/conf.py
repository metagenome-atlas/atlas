# Configuration file for the Sphinx documentation builder.

# -- Project information

project = "Metagenome-atlas"
copyright = "2021, Silas Kieser"
author = "Joe Brown and Silas Kieser"
# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '2.8'
# The full version, including alpha/beta/rc tags.
release = '2.8.0'


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ["build", "Thumbs.db", ".DS_Store", "old"]

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
