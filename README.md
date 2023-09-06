<!---

**⚠ Warning**

```This code is in a pre-release state. It may not be fully functional and breaking changes can occur without notice.```

--->

# Gurobi Model Analyzer

Gurobi Model Analyzer is an
[open-source](https://gurobi-modelanalyzer.readthedocs.io/en/stable/license.html) python package that provides
detailed analysis of model solutions and model characteristics.
Version 1.0 consists of a results_analyzer module that calculates
explanations of ill-conditioned basis matrices.

Note: At the moment this repository is public, but final docs configuration
and packages is not yet finished.   We recommend that you hold off on using it
until version 1.0 is announced as available.

# Documentation

The latest user manual is available on
[readthedocs](https://gurobi-optimization-gurobi-modelanalyzer.readthedocs-hosted.com/en/latest/).


# Contact us

For questions related to using Gurobi Model Analyzer, please use
[Gurobi's Forum](https://support.gurobi.com/hc/en-us/community/topics/10373864542609-GitHub-Projects).

For reporting bugs, issues, and feature requests please open an issue.

If you encounter issues with Gurobi or gurobipy please contact
[Gurobi Support](https://support.gurobi.com/hc/en-us).


# Installation

## Dependencies

- Python >= 3.9
- [`numpy`](https://pypi.org/project/numpy/)  >= 1.21.5    (although earlier
  versions compatible with python 3.7 will probably work).


## Pip installation

The easiest way to install gurobi-modelanalyzer is using pip in a
virtual environment:

```
(.venv) pip install gurobi-modelanalyzer
```

This will also install the numpy and gurobipy dependencies.  One of the
advanced functions makes use of matplotlib; if you haven't already installed
that and plan to use this function (matrix_bitmap), you can either install
the matplotlib package directly, or install it with the gurobi-modelanalyzer
package via "pip install gurobi-modelanalyzer matplotlib".

Please note that gurobipy is commercial software and requires a
license. When installed via pip or conda, gurobipy ships with a free
license for testing and can only solve models of limited size.


Then use the explainer functions.   Example usage

```
import gurobipy as gp
import gurobi_modelanalyzer as gma
model=gp.read("myillconditionedmodel.mps")
model.optimize()
gma.kappa_explain(model)

# row-based explanation (default)
gma.kappa_explain(model, expltype="ROWS")

# column-based explanation
gma.kappa_explain(model, expltype="COLS")

# angle-based explanation (only looks for pairs of rows or columns
# that cause ill-conditioning.
gma.angle_explain(model)
```

Use `help(gma.kappa_explain)` or `help(gma.angle_explain)` for information
on more advanced usage.


# Getting a Gurobi License
Alternatively to the bundled limited license, there are licenses that can handle models of all sizes.

As a student or staff member of an academic institution you qualify for a free, full product license.
For more information, see:

* https://www.gurobi.com/academia/academic-program-and-licenses/

For a commercial evaluation, you can
[request an evaluation license](https://www.gurobi.com/free-trial/?utm_source=internal&utm_medium=documentation&utm_campaign=fy21_pipinstall_eval_pypipointer&utm_content=c_na&utm_term=pypi).

Other useful resources to get started:
* https://www.gurobi.com/documentation/
* https://support.gurobi.com/hc/en-us/community/topics/


# Development
We value any level of experience in using Gurobi Model Analyzer and would like to encourage you to
contribute directly to this project. Please see the [Contributing Guide](CONTRIBUTING.md) for more information.

## Source code
You can clone the latest sources with the command:

```
git clone git@github.com:Gurobi/gurobi-modelanalyzer.git
```


## Testing


## Submitting a Pull Request
Before opening a Pull Request, have a look at the full
[Contributing page](CONTRIBUTING.md) to make sure your code complies with
our guidelines.
