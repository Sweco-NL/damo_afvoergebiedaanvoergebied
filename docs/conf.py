# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "DrainageUnits"
copyright = "2024, AnaisWens"
author = "AnaisWens"
release = "18/12"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

import os
import sys

sys.path.insert(
    0, os.path.abspath("../code")
)  # Pad naar de code-map relatief aan de source-folder  # Gebruik forward slashes

# Voeg extensies toe
extensions = [
    "sphinx.ext.autodoc",  # Voor het automatisch genereren van documentatie vanuit docstrings
    "sphinx.ext.napoleon",  # Voor Google en NumPy stijl docstrings
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
html_static_path = ["_static"]
