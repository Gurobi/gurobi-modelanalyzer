API Reference
#############

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
