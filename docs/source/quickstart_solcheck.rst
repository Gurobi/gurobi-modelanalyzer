Quick Start Guide
#################

When developing an optimization model, you may have solution values that you
want to test to see if everything is correct. This tool helps you do this
without writing code: it tests your solution for a given optimization model. A
full solution is not required; it can test a solution for a subset of decision
variables.


Python
******


.. code-block:: python

   import gurobipy as gp
   import gurobi_modelanalyzer as gma

   m = gp.read("examples/data/afiro.mps")

   # Provide a solution
   sol = {m.getVarByName("X01"): 78, m.getVarByName("X22"): 495}

   # Instantiate
   sc = gma.SolCheck(m)

   # Check solution
   sc.test_sol(sol)
   print(f"Solution Status: {sc.Status}")


See :doc:`usage_solcheck` section for more details and the
:doc:`apiexamples_solcheck` for some examples.

Command-Line
************

Provided a model file and a solution file (`SOL format <https://docs.gurobi.com/projects/optimizer/en/current/reference/misc/fileformats.html#sol-format>`__
or `JSON solution format <https://docs.gurobi.com/projects/optimizer/en/current/reference/misc/fileformats.html#json-solution-format>`__)
we can also use the ``gurobi_solcheck`` command-line tool:

::

   gurobi_solcheck --model examples/data/afiro.mps --sol afiro.json


See :doc:`usage_solcheck` section for more details and the
:doc:`cli_examples_solcheck` for some examples.
