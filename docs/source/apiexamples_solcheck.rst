API Examples
############


Suboptimal solutions
====================

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
===========================

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
