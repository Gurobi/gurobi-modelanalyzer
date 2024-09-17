# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


project = "Gurobi Model Analyzer"
copyright = "2024, Gurobi Optimization"
author = "Gurobi Optimization"
html_title = "Gurobi Model Analyzer"

html_theme = "gurobi_sphinxtheme"
html_favicon = "https://www.gurobi.com/favicon.ico"

# -- Warning banner while in beta -------------------------------------------
# rst_prolog = """.. warning::
#    This code is in a pre-release state. It may not be fully functional and breaking changes
#    can occur without notice.
# """

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
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

html_static_path = ["_static"]
