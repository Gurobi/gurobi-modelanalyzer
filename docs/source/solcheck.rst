Solution checker
================

Test a known solution in an optimization model using Gurobi Optimizer

Background
----------

When developing an optimization model, you may have solution values that
you want to test to see if everything is correct. This tool helps you do
this without writing code: it tests your solution for a given
optimization model. A full solution is not required; it can test a
solution for a subset of decision variables.

Examples
--------

Solution is feasible but suboptimal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   gurobi_solcheck --model examples/data/afiro.mps --sol afiro.json

This gives the following output:

::

   Objectives:
   Fixed:      -459.6919
   Optimal:    -464.7531
   Difference: -5.0613

So the solution file afiro.sol is -5.0613 worse than an optimal
solution. We can get an optimal solution by adding the ``--result``
flag:

::

   gurobi_solcheck --model examples/data/afiro.mps --sol afiro.sol --result afirofix

This writes the solution file afirofix.sol, which we can compare with
the solution file afiro.sol.

Adjusting constraints to make an infeasible solution feasible
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   gurobi_solcheck --model examples/data/misc07.mps --sol misc07.sol

This gives the following output:


.. code-block:: none

   Solution is infeasible for feasibility tolerance of 1e-06
   Fixed values are 1.5 from a feasible solution

By default, this finds the minimum adjustment to the right hand sides of
the constraints to make the solution feasible. We can understand this
better by writing a result file:

::

   gurobi_solcheck --model examples/data/misc07.mps --sol misc07.sol --result misc07fix

This produces the violation file misc07fix.vio, which has the following
nonzero values:

::

   ROW001 -0.5
   ROW074 1.0

This means the solution would be feasible if the right hand side of
constraint ROW001 is reduced by 0.5 and the right hand side of
constraint of constraint ROW074 is increased by 1.0; the total
adjustment is 1.5 as reported in the log.

Adjusting solution values to make an infeasible solution feasible
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can adjust the solution values to make it feasible by setting the
flag ``--infmethod V``:

::

   gurobi_solcheck --model examples/data/misc07.mps --sol misc07.sol --infmethod V --result misc07fix

This produces the solution file misc07fix.sol. Comparing misc07.sol and
misc07fix.sol, we see that COL260 is 2400.5 in misc07.sol and 2810 in
misc07fix.sol, meaning that COL260 can be increased by 409.5 to make the
solution feasible.

Explaining an infeasible solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If we set the flag ``--infmethod I`` and specify a result file, we can
get an explanation for why the solution file is infeasible:

::

   gurobi_solcheck --model examples/data/misc07.mps --sol misc07.sol --infmethod I --result misc07fix

This produces the smaller model file misc07fix.ilp, which contains one
constraint and one fixed bound (output condensed for readability):

::

   Minimize

   Subject To
    ROW001: - 100 COL004 - 100 COL005 - 100 COL006 - 110 COL007 - 110 COL008
      - 110 COL009 - 495 COL010 - 495 COL011 - 495 COL012 - 445 COL013
      - 445 COL014 - 445 COL015 - 300 COL016 - 300 COL017 - 300 COL018
        [...]
      - 75 COL248 - 75 COL249 - 75 COL250 - 75 COL251 - 75 COL252 - 75 COL253
      - 75 COL254 - 75 COL255 - 75 COL256 - 75 COL257 - 75 COL258 - 75 COL259
      + COL260 = 0
   Bounds
    [...]
    COL260 = 2400.5
   Generals
    COL004 COL005 COL006 COL007 COL008 COL009 COL010 COL011 COL012 COL013
    [...]
    COL258 COL259
   End

Here, all variables besides COL260 are integers, so constraint ROW001
cannot be feasible for COL260 = 2400.5.

How it works
------------

When testing a solution for an optimization model, there are several
possible cases:


.. raw:: html

   <ol>
   <li>The solution is feasible but possibly suboptimal</li>
   <li>The solution is infeasible and:<br/>
   <span>&nbsp;&nbsp; a. You want an explanation why it's infeasible</span><br/>
   <span>&nbsp;&nbsp; b. You want to repair the solution</span><br/>
   <span>&nbsp;&nbsp; c. You want to repair the model</span><br/>
   </li>
   </ol>


This tool can diagnose all of these cases. Note that the difference
between cases 2b and 2c is a policy decision: in case of infeasibility,
you need to determine whether you want to repair the solution or the
model.

Case 1: The solution is feasible but possibly suboptimal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, the solution checker verifies that the given solution satisfies
the constraints. If so, it will also compute an optimal solution, to
show how far the given solution is from optimal. If you write a .sol
solution file, you can see the difference from the given solution and an
optimal solution.

Case 2a: The solution is infeasible and you want an explanation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the solution checker verifies that the given solution does not
satisfy the model constraints, you can use it to find a set of
constraints that conflicts with the given solution. This is done by
setting *infmethod* to I, which calls IIS to find a set of constraints
that conflicts with the given solution. Note that if there are multiple
independent sets of constraints that conflict, this will only find one
set of conflicts; to check for others, remove the infeasible constraints
from the model and try again.

Case 2b: The solution is infeasible and you want to repair the solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the solution checker verifies that the given solution does not
satisfy the model constraints, you can use it to modify the solution to
make it feasible. This is done by setting *infmethod* to V, which calls
feasRelax to find the smallest changes to the solution values that make
the constraints feasible.

Case 2c: The solution is infeasible and you want to repair the model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the solution checker verifies that the given solution does not
satisfy the model constraints, you can use it to modify the solution to
make it feasible. This is done by setting *infmethod* to C, which calls
feasRelax to find the smallest changes to the constraints that make the
solution feasible.

Usage
-----

Interactive mode
~~~~~~~~~~~~~~~~

If you run gurobi_solcheck without specifying a model file, it will run
in interactive mode. In interactive mode, it will prompt for a model
file and a solution file (or zero values). It will also prompt whether
to diagnose the solution, and for an optional result file.

Example:

::

   Entering interactive mode
   ? Select a model file /Library/gurobi1103/macos_universal2/examples/data/afiro.mps
   Read MPS format model from file /Library/gurobi1103/macos_universal2/examples/data/afiro.mps
   Reading time = 0.00 seconds
   AFIRO: 27 rows, 32 columns, 83 nonzeros
   ? Use a solution file? Yes
   ? Select a solution file ../afiro.sol
   ? Ignore zeros in solution file? Yes

   [...]

   Solution is feasible for feasibility tolerance of 1e-06

   ? Diagnose solution? (Y/n)

Command-Line
~~~~~~~~~~~~

You can run gurobi_solcheck via the command-line using the following
syntax:

::

   gurobi_solcheck --model [modelfile] {--sol [solutionfile]} {--testonly} {--infmethod [C,V,I]} {--result [resultfile]}

This reads the model file (.mps, .lp, etc.) and tests it with the
solution file (`SOL format <https://docs.gurobi.com/projects/optimizer/en/current/reference/misc/fileformats.html#sol-format>`__
or `JSON solution format <https://docs.gurobi.com/projects/optimizer/en/current/reference/misc/fileformats.html#json-solution-format>`__).
The solution file can be optionally compressed with gzip (.sol.gz). If
you omit the solution file, it will test the solution of all zeros.

To simply test if the solution is feasible and skip diagnostics, add the
flag ``--testonly``.

If a solution is infeasible, use ``--infmethod`` to select the method to
diagnose the infeasible solution:

-  C: repairs by adjusting the right-hand-side of constraints (default)
-  V: repairs by finding feasible variable values that are close to the
   solution values
-  I: computes the `Irreducible Inconsistent Subsystem <https://docs.gurobi.com/projects/optimizer/en/current/reference/python/model.html#Model.computeIIS>`__
   to explain an infeasible model

If you want to save the result, add the ``--result`` flag to specify the
result file, which will be saved in .sol format for a solution, .vio
format for constraint violations and .ilp for an IIS. The .vio format is
similar to the .sol format, except that it uses constraint names.

The exit status code will be the `solution status <https://docs.gurobi.com/projects/optimizer/en/current/reference/attributes/model.html#status>`__
of the specified solution. If you skip diagnostics and the solution is
feasible, the status will be GRB.SUBOPTIMAL, even if the solution is
optimal.

Command-Line examples
~~~~~~~~~~~~~~~~~~~~~

::

   gurobi_solcheck --model examples/data/afiro.mps --sol afiro.sol
   gurobi_solcheck --model examples/data/misc07.mps --sol misc07.sol --result misc07fix

Gurobi Parameters
~~~~~~~~~~~~~~~~~

If you want to set gurobi parameters with the command-line tool, create
a file called ``gurobi.env`` in your working directory; for details, see
the section "Using a gurobi.env file" in `this Gurobi knowledge base article about setting parameters <https://support.gurobi.com/hc/en-us/articles/14126085438481-How-do-I-set-parameters-in-Gurobi>`__.

API Reference
-------------

Everything is done using the ``SolCheck`` class:


.. py:class:: SolCheck

  .. py:method:: SolCheck(model)

     Initializes a SolCheck using a ``model`` a ``gurobipy.Model`` object.


  .. py:method:: SolCheck.test_sol(sol)

     Test the solution values ``sol``, a Python dictionary where the keys are
     gurobipy.Var objects and the values are the solution values. This only
     tests if the solution values are feasible or not; you must call
     additional methods to diagnose the solution values.


  .. py:method:: SolCheck.inf_explain()

     Computes the `Irreducible Inconsistent Subsystem <https://docs.gurobi.com/projects/optimizer/en/current/reference/python/model.html#Model.computeIIS>`__
     to explain an infeasible solution.

  .. py:method:: SolCheck.inf_repair(repairMethod='C', makeCopy=False)

     Repairs an infeasible solution.

     :param repairMethod: String to set the method to use to repair the
       infeasibility. If it is ``"C"`` (default), it repairs by adjusting the
       right-hand-side values of constraints; if it  is ``"V"``, it repairs by
       adjusting the solution values.

     :param makeCopy: Bool to make a fresh copy of the model object; if it is
       ``False`` (defult), then the original model object will be replaced by the
       relaxed copy.

     For infeasible models where repairMethod is ``"C"``, the ``Constr`` objects
     will have an additional floating point attribute ``_Violation`` that measures
     how much that constraint is violated; ``_Violation`` may be positive or
     negative.


  .. py:method:: SolCheck.optimize()

     Optimizes the original model, starting from the test solution. For a solution
     that is feasible, this can determine how far that solution may be from
     optimal.


  .. py:method:: SolCheck.write_result(fn)

     If you call any of the explanation methods (``SolCheck.inf_explain()``,
     ``SolCheck.inf_repair()`` or ``SolCheck.optimize()``), this will write a
     result file; the type of result will depend on the solution status and the
     type of explanation.


Examples
~~~~~~~~


Suboptimal solutions
^^^^^^^^^^^^^^^^^^^^

Given a solution, we can use the :py:meth:`SolCheck.test_sol` function to test
it.

.. code-block:: python


    import gurobipy as gp
    import gurobi_modelanalyzer as gma

    m = gp.read("afiro.mps")

    sol = {m.getVarByName("X01"): 78, m.getVarByName("X22"): 495}
    sc = gma.SolCheck(m)

    sc.test_sol(sol)
    print(f"Solution Status: {sc.Status}")


Will print:

.. code-block:: none

   Read MPS format model from file afiro.mps
   Reading time = 0.00 seconds
   AFIRO: 27 rows, 32 columns, 83 nonzeros
   Gurobi Optimizer version 11.0.3 build v11.0.3rc0

   Thread count: 8 physical cores, 8 logical processors, using up to 8 threads

   Optimize a model with 27 rows, 32 columns and 83 nonzeros
   Model fingerprint: 0x540c3b7f
   Coefficient statistics:
     Matrix range     [1e-01, 2e+00]
     Objective range  [3e-01, 1e+01]
     Bounds range     [8e+01, 5e+02]
     RHS range        [4e+01, 5e+02]
   Presolve removed 19 rows and 22 columns
   Presolve time: 0.03s
   Presolved: 8 rows, 10 columns, 27 nonzeros

   Iteration    Objective       Primal Inf.    Dual Inf.      Time
          0   -4.5969189e+02   2.146875e+00   0.000000e+00      0s
          3   -4.5969189e+02   0.000000e+00   0.000000e+00      0s

   Solved in 3 iterations and 0.05 seconds (0.00 work units)
   Optimal objective -4.596918857e+02

   Solution is feasible for feasibility tolerance of 1e-06

   Solution Status: 13



We can check this solution against the optimal one by calling
:py:meth:`SolCheck.optimize`.

.. code-block:: python

   sc.optimize()
   for v in sol.keys():
       print(f"{v.VarName}: Fixed value: {sol[v]}, Computed value: {v.X}")

Produces

.. code-block:: none

   Comparing quality with original solution

   Gurobi Optimizer version 11.0.3 build v11.0.3rc0

   Thread count: 8 physical cores, 8 logical processors, using up to 8 threads

   Optimize a model with 27 rows, 32 columns and 83 nonzeros
   Coefficient statistics:
     Matrix range     [1e-01, 2e+00]
     Objective range  [3e-01, 1e+01]
     Bounds range     [0e+00, 0e+00]
     RHS range        [4e+01, 5e+02]

   Iteration    Objective       Primal Inf.    Dual Inf.      Time
          0   -3.1277714e+30   1.240950e+31   3.127771e+00      0s
          4   -4.6475314e+02   0.000000e+00   0.000000e+00      0s

   Solved in 4 iterations and 0.00 seconds (0.00 work units)
   Optimal objective -4.647531429e+02

   Objectives:
   Fixed:      -459.6919
   Optimal:    -464.7531
   Difference: -5.0613

   X01: Fixed value: 78, Computed value: 80.0
   X22: Fixed value: 495, Computed value: 500.0


We can see that the solution we provided is worse than the optimal solution by
-5.0613 in total, and the difference in the solution values that we provided.

Test an infeasible solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    m = gp.read("misc07.mps")
    sol = {m.getVarByName("COL260"): 2400.5}
    sc = gma.SolCheck(m)

    sc.test_sol(sol)

    print(f"Solution Status: {sc.Status}")

Will print:

.. code-block:: none

   Read MPS format model from file misc07.mps
   Reading time = 0.00 seconds
   MISC07: 212 rows, 260 columns, 8619 nonzeros
   Gurobi Optimizer version 11.0.3 build v11.0.3rc0

   Thread count: 8 physical cores, 8 logical processors, using up to 8 threads

   Optimize a model with 212 rows, 260 columns and 8619 nonzeros
   Model fingerprint: 0xd79ad074
   Variable types: 1 continuous, 259 integer (0 binary)
   Coefficient statistics:
     Matrix range     [1e+00, 7e+02]
     Objective range  [1e+00, 1e+00]
     Bounds range     [1e+00, 2e+03]
     RHS range        [1e+00, 3e+02]
   Presolve removed 0 rows and 7 columns
   Presolve time: 0.00s

   Explored 0 nodes (0 simplex iterations) in 1.72 seconds (0.00 work units)
   Thread count was 1 (of 8 available processors)

   Solution count 0

   Model is infeasible
   Best objective -, best bound -, gap -

   Solution is infeasible for feasibility tolerance of 1e-06

   Solution Status: 3

Here we can see that the solution we provided makes the problem infeasible. We
can use the :py:meth:`SolCheck.inf_repair` function to repair the infeasibility.

.. code-block:: python

   sc.inf_repair()
   for c in m.getConstrs():
       if abs(c._Violation) > 0.0001:
           print(f"{c.ConstrName}: RHS: {c.RHS}, Violation: {c._Violation}")

We get:

.. code-block:: none

   Relaxing to find smallest violation from fixed solution

   Gurobi Optimizer version 11.0.3 build v11.0.3rc0

   Thread count: 8 physical cores, 8 logical processors, using up to 8 threads

   Optimize a model with 212 rows, 507 columns and 8866 nonzeros
   Model fingerprint: 0x396303c6
   Variable types: 248 continuous, 259 integer (0 binary)
   Coefficient statistics:
     Matrix range     [1e+00, 7e+02]
     Objective range  [1e+00, 1e+00]
     Bounds range     [1e+00, 2e+03]
     RHS range        [1e+00, 3e+02]
   Found heuristic solution: objective 2534.5000000
   Presolve removed 0 rows and 7 columns
   Presolve time: 0.01s
   Presolved: 212 rows, 500 columns, 8823 nonzeros
   Variable types: 70 continuous, 430 integer (380 binary)

   Root relaxation: objective 0.000000e+00, 124 iterations, 0.00 seconds (0.00 work units)

       Nodes    |    Current Node    |     Objective Bounds      |     Work
    Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

        0     0    0.00000    0   22 2534.50000    0.00000   100%     -    0s
   H    0     0                     207.5000000    0.00000   100%     -    0s
   H    0     0                     173.5000000    0.00000   100%     -    0s
   H    0     0                     110.5000000    0.00000   100%     -    0s
   H    0     0                      97.5000000    0.00000   100%     -    0s
   H    0     0                      61.5000000    0.00000   100%     -    0s
   H    0     0                      12.5000000    0.00000   100%     -    0s
        0     0    0.50000    0   25   12.50000    0.50000  96.0%     -    0s
   H    0     0                       7.5000000    0.50000  93.3%     -    0s
   H    0     0                       6.5000000    0.50000  92.3%     -    0s
        0     0    0.50000    0   34    6.50000    0.50000  92.3%     -    0s
        0     0    0.50000    0   31    6.50000    0.50000  92.3%     -    0s
        0     0    0.50000    0   30    6.50000    0.50000  92.3%     -    0s
        0     0    0.50000    0   27    6.50000    0.50000  92.3%     -    0s
        0     0    0.50000    0   23    6.50000    0.50000  92.3%     -    0s
        0     0    0.50000    0   27    6.50000    0.50000  92.3%     -    0s
        0     0    0.50000    0   22    6.50000    0.50000  92.3%     -    0s
   H    0     0                       5.5000000    0.50000  90.9%     -    0s
        0     2    0.50000    0   22    5.50000    0.50000  90.9%     -    0s
   H   80    88                       4.5000000    0.50000  88.9%  41.2    0s
   H  132   209                       3.5000000    0.50000  85.7%  37.7    0s
   *  706   534              39       2.5000000    0.50000  80.0%  24.9    0s
   H 1487   801                       2.5000000    0.50000  80.0%  25.0    1s
   H 1490   801                       2.5000000    0.50000  80.0%  25.0    1s
     5201  1713    0.58621   24   19    2.50000    0.50000  80.0%  24.9    5s
   * 6181   908              26       1.5000000    0.50000  66.7%  24.6    5s

   Cutting planes:
     Gomory: 4
     MIR: 8
     Flow cover: 65

   Explored 9399 nodes (234850 simplex iterations) in 6.19 seconds (6.17 work units)
   Thread count was 8 (of 8 available processors)

   Solution count 10: 1.5 2.5 2.5 ... 61.5

   Optimal solution found (tolerance 1.00e-04)
   Best objective 1.500000000000e+00, best bound 1.500000000000e+00, gap 0.0000%

   Fixed values are 1.5 from a feasible solution

   ROW001: RHS: 0.0, Violation: -0.5
   ROW074: RHS: 1.0, Violation: 1.0


From this we can see that we would have to relax constraints ``ROW001`` and
``ROW074`` by -0.5 and 1.0 to make the problem feasible.

Tips and Tricks
---------------

Handling slow runtime
~~~~~~~~~~~~~~~~~~~~~

By default, this solution checker has no time limit. In some cases, it
may take lots of time and memory to test a solution. In the worst case,
this could be comparable to solving the model itself from scratch.

If runtime is an issue, try the following:

1. Test a larger subset of variables. Fixing more variables leaves a
   smaller, easier test problem

2. Test zero values. By default, the command-line version ignores zero
   values, but you can force it to test zero values by adding the
   ``--usezeros`` flag.

3. Set a runtime limit. While this may not produce the ideal explanation
   or repair, it may still give useful information. For the command-line
   version, set the TimeLimit parameter using the instructions above on
   `Gurobi Parameters <#Gurobi-Parameters>`__.

About tolerances for infeasible models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tool uses the Gurobi solver to test solutions for a model. In some
cases, a solution may be on the border between feasible and infeasible.
This depends on the value of FeasibilityTol, the feasibility tolerance
parameter. If the solution is infeasible, the solution checker will
report the current value of FeasibilityTol. You may want to increase
FeasibilityTol to test if the solution becomes feasible. To do this with
the command-line, follow the instructions above on `Gurobi
Parameters <#Gurobi-Parameters>`__.

For an integer program, a near-integer solution value like 0.9999 may be
rejected as infeasible. We recommend using precise solution values like
1.0.
