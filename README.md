# Model analyzer

The ill-conditioning explainer, packaged! Install using:

```
pip install git+ssh://git@github.compute.gurobi.com/bowly/model_analyzer.git
```

Then use the explainer functions:

```
from model_analyzer import kappa_explain, BYROWS, BYCOLS

# row-based explanation
kappa_explain(model, expltype=BYROWS)

# column-based explanation
kappa_explain(model, expltype=BYCOLS)
```

I have only exposed the `kappa_explain` function for now (ill-conditioning checker)
as it's not yet clear what else will be in the public API.
