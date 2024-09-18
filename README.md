# Gurobi Model Analyzer

Gurobi Model Analyzer is an
[open-source](https://gurobi-modelanalyzer.readthedocs.io/en/stable/license.html) python package that provides
detailed analysis of model solutions and model characteristics.
It consists of a results_analyzer module that calculates
explanations of ill-conditioned basis matrices and a solcheck module that analysizes a given solution.


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
- [questionary](https://pypi.org/project/questionary/) for the (optional) interactive version

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


## Example usage
### Using the explainer functions

```python
import gurobipy as gp
import gurobi_modelanalyzer as gma

model = gp.read("myillconditionedmodel.mps")
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

### Using the solution checker

Testing a suboptimal solution

```python
import gurobipy as gp
import gurobi_modelanalyzer as gma

m = gp.read("examples/data/afiro.mps")

sol = {m.getVarByName("X01"): 78, m.getVarByName("X22"): 495}
sc = gma.SolCheck(m)

sc.test_sol(sol)
print(f"Solution Status: {sc.Status}")
sc.optimize()
for v in sol.keys():
    print(f"{v.VarName}: Fixed value: {sol[v]}, Computed value: {v.X}")
```

Testing an infeasible solution

```python
m = gp.read("examples/data/misc07.mps")

sol = {m.getVarByName("COL260"): 2400.5}
sc = gma.sol_check(m)

sc.test_sol(sol)

print(f"Solution Status: {sc.Status}")
sc.inf_repair()
for c in m.getConstrs():
    if abs(c._Violation) > 0.0001:
        print(f"{c.ConstrName}: RHS: {c.RHS}, Violation: {c._Violation}")
```


# Getting a Gurobi License
Alternatively to the bundled limited license, there are licenses that can handle models of all sizes.

As a student or staff member of an academic institution, you qualify for a free, full-product license.
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


## Submitting a Pull Request
Before opening a Pull Request, have a look at the full
[Contributing page](CONTRIBUTING.md) to make sure your code complies with
our guidelines.
