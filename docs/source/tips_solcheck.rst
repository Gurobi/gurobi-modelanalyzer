Tips and Tricks
###############

Handling slow runtime
*********************

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
**************************************

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

Gurobi Parameters
*****************

If you want to set gurobi parameters with the command-line tool, create
a file called ``gurobi.env`` in your working directory; for details, see
the section "Using a gurobi.env file" in `this Gurobi knowledge base article about setting parameters <https://support.gurobi.com/hc/en-us/articles/14126085438481-How-do-I-set-parameters-in-Gurobi>`__.
