.. _ScalingCLILabel:

Command-Line Interface
######################

The ``gurobi_cls`` command-line tool scales a Gurobi model and solves it in
one step, mimicking the interface of Gurobi's own ``gurobi_cl`` tool.


Basic Usage
***********

.. code-block:: none

   gurobi_cls [scaling options] [Param=Value ...] model

The only required argument is the path to the model file. Any Gurobi solver
parameters can be passed as ``Param=Value`` positional arguments before the
model path, exactly as with ``gurobi_cl``:

.. code-block:: none

   gurobi_cls --method equilibration TimeLimit=60 Presolve=2 model.mps


Output
******

``gurobi_cls`` produces three outputs:

1. **Console**: the scaling log and, after solving, solution quality statistics
   for both the scaled and the original (unscaled) variable space.

2. **Log file** (default: ``gurobi.log``): a unified file containing the
   scaling log, the Gurobi solve log, and the solution quality block. Override
   the name with ``LogFile=path``.

3. **Scaling file** (default: ``<model_stem>.scl``): the computed scaling
   factors in ``.scl`` format, written automatically after scaling completes.
   All factors are written with ``lock_flag=0`` so that re-importing the file
   reproduces the same scaling exactly.


Scaling Options
***************

.. option:: --method {equilibration,geometric_mean,arithmetic_mean}, -m {equilibration,geometric_mean,arithmetic_mean}

   Scaling algorithm to apply. Default: ``equilibration``.

   .. list-table::
      :header-rows: 1
      :widths: 25 11 11 11 11 11

      * - Method
        - (MI)LP
        - (MI)QP
        - (MI)QCP
        - (MI)QCQP
        - (MI)NLP
      * - ``equilibration``
        - ✓
        - ✓
        - ✓
        - ✓
        - —
      * - ``geometric_mean``
        - ✓
        - —
        - ✓
        - —
        - —
      * - ``arithmetic_mean``
        - ✓
        - —
        - ✓
        - —
        - —

   (MI)QCP denotes quadratic constraints with a linear objective; (MI)QCQP
   denotes quadratic constraints and a quadratic objective; see
   :ref:`ScalingAdvUsageLabel` for details. (MI)NLP models are not currently
   supported.

.. option:: --scale-passes N

   Maximum number of scaling iterations. Default: 1.

.. option:: --scale-conv-tol TOL

   Convergence tolerance. Scaling stops early when the maximum deviation of
   the scaling factors from 1 falls below this value. Default: 1e-4.

.. option:: --scaling-lb LB

   Lower bound for scaling factors. Default: 1e-8.

.. option:: --scaling-ub UB

   Upper bound for scaling factors. Default: 1e8.

.. option:: --value-threshold THRESH

   Coefficients with absolute value below this threshold are treated as zero
   after scaling. Default: 1e-13.

.. option:: --scaling-time-limit SEC

   Wall-clock time limit in seconds for the scaling iterations. If reached,
   the best scaling found so far is used and the algorithm stops. Default: no
   limit.

.. option:: --scaling-file PATH

   Path to a ``.scl`` scaling input file. When provided, the factors in the
   file are pre-applied to the model before the scaling algorithm runs
   (``init_scaling=2`` warmstart mode). See :ref:`ScalingFilesLabel` for the
   file format.

   This flag is the CLI counterpart of
   :py:func:`~gurobi_modelanalyzer.read_scaling_file`.

.. option:: --no-console-log

   Suppress the scaling log on the console. The log is still written to the
   log file.


Exit Codes
**********

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Code
     - Meaning
   * - 0
     - Optimal or suboptimal solution found.
   * - 1
     - Error reading the model or scaling file.
   * - 2
     - Scaled model is infeasible.
   * - 3
     - Scaled model is unbounded.
   * - 4
     - Solver finished with a non-optimal status.


Examples
********

Scale and solve with default settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: none

   gurobi_cls model.mps

This scales ``model.mps`` with the equilibration method (1 pass), solves the
scaled model, writes the solution quality to the console, and saves the log to
``gurobi.log`` and the scaling factors to ``model.scl``.


Pass Gurobi solver parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: none

   gurobi_cls TimeLimit=120 Method=2 LogFile=run.log model.mps

Parameters are forwarded directly to the Gurobi solver. ``LogFile`` controls
the name of the unified output log.


Use a different scaling method and tighter convergence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: none

   gurobi_cls --method geometric_mean --scale-passes 10 --scale-conv-tol 1e-6 model.mps


Re-apply previously saved scaling factors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After a first run that produced ``model.scl``, the same scaling can be
reproduced exactly on a later run:

.. code-block:: none

   gurobi_cls --scaling-file model.scl model.mps

Because the factors in the file carry ``lock_flag=0``, the algorithm treats
them as fixed and no further modification occurs.


Use scaling factors as a warmstart for a related model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scaling factors from one model can serve as a starting point for a related
model. Edit the ``.scl`` file to change ``lock_flag`` from ``0`` to ``1`` on
the entries that should be refined, then run:

.. code-block:: none

   gurobi_cls --scaling-file previous.scl related_model.mps

The algorithm will start from those factors and continue iterating.


Suppress console output and redirect the log
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: none

   gurobi_cls --no-console-log LogFile=nightly.log model.mps


Python API / CLI Round-Trip
***************************

Scaling factors can flow freely between the Python API and the command-line
tool via ``.scl`` files:

.. code-block:: python

   # Compute scaling in Python, then run the CLI with the same factors
   import gurobipy as gp
   import gurobi_modelanalyzer as gma

   m = gp.read("model.mps")
   m_scaled = gma.scale_model(m, method="equilibration")
   m_scaled.write_scaling("model.scl")

.. code-block:: none

   # Reproduce the same scaling via the CLI
   gurobi_cls --scaling-file model.scl model.mps

Conversely, factors produced by ``gurobi_cls`` (written automatically to
``model.scl``) can be loaded in Python:

.. code-block:: python

   from gurobi_modelanalyzer.scaling import scale_model, read_scaling_file

   m = gp.read("model.mps")
   read_scaling_file("model.scl", m)
   m_scaled = scale_model(m, method="equilibration", init_scaling=2)
