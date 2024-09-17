How it works
############

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
********************************************************

First, the solution checker verifies that the given solution satisfies
the constraints. If so, it will also compute an optimal solution, to
show how far the given solution is from optimal. If you write a .sol
solution file, you can see the difference from the given solution and an
optimal solution.

Case 2a: The solution is infeasible and you want an explanation
***************************************************************

If the solution checker verifies that the given solution does not
satisfy the model constraints, you can use it to find a set of
constraints that conflicts with the given solution. This is done by
setting *infmethod* to I, which calls IIS to find a set of constraints
that conflicts with the given solution. Note that if there are multiple
independent sets of constraints that conflict, this will only find one
set of conflicts; to check for others, remove the infeasible constraints
from the model and try again.

Case 2b: The solution is infeasible and you want to repair the solution
***********************************************************************

If the solution checker verifies that the given solution does not
satisfy the model constraints, you can use it to modify the solution to
make it feasible. This is done by setting *infmethod* to V, which calls
feasRelax to find the smallest changes to the solution values that make
the constraints feasible.

Case 2c: The solution is infeasible and you want to repair the model
********************************************************************

If the solution checker verifies that the given solution does not
satisfy the model constraints, you can use it to modify the solution to
make it feasible. This is done by setting *infmethod* to C, which calls
feasRelax to find the smallest changes to the constraints that make the
solution feasible.
