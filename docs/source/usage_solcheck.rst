Usage
#####

Interactive mode
****************

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
************

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


Python API
**********

All the functionality can also be accessed using the :py:class:`SolCheck`
object.


Examples
########

.. toctree::
   :maxdepth: 2

   cli_examples_solcheck
   apiexamples_solcheck
