# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Gurobi Model Analyzer"
copyright = "2023, Gurobi Optimization, LLC. All Rights Reserved."
html_logo = "_static/gurobi-logo-title.png"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.extlinks",
    "sphinxcontrib.bibtex",
]

bibtex_bibfiles = ["bib_illcond.bib"]

extlinks_detect_hardcoded_links = True
extlinks = {
    "issue": ("https://github.com/Gurobi/gurobi-machinelearning/issues/%s", "issue %s"),
    "gurobipy": (
        "https://www.gurobi.com/documentation/current/refman/py_%s.html",
        "gurobipy %s",
    ),
    "pypi": ("https://pypi.org/project/%s/", "%s"),
}


templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
