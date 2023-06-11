# Model analyzer

Gurobi Model Analyzer

Gurobi Model Analyzer is an open sources python package that provides
detailed analysis of model solutions and model characteristics.
Version 1.0 consists of a results_analyzer module that calculates
explanations of ill conditioned basis matrices.


Documentation

The latest user manual is available on readthedocs


Contact us

For questions related to using Gurobi Machine Learning please use
Gurobi's Forum.

For reporting bugs, issues and feature request please open an issue.

If you encounter issues with Gurobi or gurobipy please contact Gurobi Support.


Installation


Dependencies

numpy >= 1.21.5    (although earlier versions compatible with python 3.7 will
                    probably work).


Pip installation

```

The easiest way to install gurobi-machinelearning is using pip in a
virtual environment:

<this needs to updated>
(.venv) pip install git+ssh://git@github.compute.gurobi.com/bowly/model_analyzer.git



This will also install the numpy and gurobipy dependencies.

Please note that gurobipy is commercial software and requires a
license. When installed via pip or conda, gurobipy ships with a free
license which is only for testing and can only solve models of limited
size.

```

Then use the explainer functions.   Example usage

```
import gurobipy as gp
import model_analyzer.results_analyzer as ra
model=gp.read("myillconditionedmodel.mps")
model.optimize()
ra.kappa_explain(model)

# row-based explanation (default)
ra.kappa_explain(model, expltype="ROWS")

# column-based explanation
ra.kappa_explain(model, expltype="COLS")

# angle-based explanation (only looks for pairs of rows or columnns
# that cause ill conditioning.

ra.angle_explain(model)
```
Use help(ra.kappa_explain) or help(ra.angle_explain) for information
on more advanced usage.