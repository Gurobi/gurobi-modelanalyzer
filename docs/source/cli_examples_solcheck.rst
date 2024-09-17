CLI Examples
############

Solution is feasible but suboptimal
***********************************

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
*************************************************************

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
*****************************************************************


We can adjust the solution values to make it feasible by setting the
flag ``--infmethod V``:

::

   gurobi_solcheck --model examples/data/misc07.mps --sol misc07.sol --infmethod V --result misc07fix

This produces the solution file misc07fix.sol. Comparing misc07.sol and
misc07fix.sol, we see that COL260 is 2400.5 in misc07.sol and 2810 in
misc07fix.sol, meaning that COL260 can be increased by 409.5 to make the
solution feasible.

Explaining an infeasible solution
*********************************

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
