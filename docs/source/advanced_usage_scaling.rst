.. _ScalingAdvUsageLabel:

Advanced Usage Guide
####################

Introduction
************

The :ref:`ScalingQSGuideLabel` covers the basic workflow for scaling a model
and recovering the unscaled solution. This section describes the more advanced
capabilities of the Model Scaling module: choosing between scaling methods,
controlling convergence, opting out selected variables or constraints,
providing user-defined initial scaling factors, accessing the scaling matrices
directly, and computing the unscaled objective value.


Choosing a Scaling Method
*************************

Three iterative scaling methods are available. The table below summarises
which model types each supports.

.. note::

   This module uses a finer classification of model types than Gurobi's own
   documentation. Gurobi refers to models with quadratic constraints generically
   as "QCP", regardless of whether a quadratic objective is also present. Here
   we distinguish:

   * **(MI)QCP**: quadratic constraints with a **linear** objective.
   * **(MI)QCQP**: quadratic constraints **and** a quadratic objective.

   This distinction matters because support for the three scaling methods
   depends on whether a quadratic objective is present, not on whether
   quadratic constraints are present.

+--------------------+----------+----------+----------+----------+----------+
| Method             | (MI)LP   | (MI)QP   | (MI)QCP  | (MI)QCQP | (MI)NLP  |
+====================+==========+==========+==========+==========+==========+
| ``equilibration``  |  ✓       |  ✓       |  ✓       |  ✓       |  —       |
+--------------------+----------+----------+----------+----------+----------+
| ``geometric_mean`` |  ✓       |  —       |  ✓       |  —       |  —       |
+--------------------+----------+----------+----------+----------+----------+
| ``arithmetic_mean``|  ✓       |  —       |  ✓       |  —       |  —       |
+--------------------+----------+----------+----------+----------+----------+

(MI)NLP models are not currently supported by any scaling method.

**Equilibration** is the recommended default. It iteratively scales rows and
columns to bring the magnitudes of the nonzero coefficients to a similar range,
following the approach described in :cite:`Elble2011`. It is the only method
that supports models with a quadratic objective ((MI)QP and (MI)QCQP). For
(MI)LP and (MI)QCP models (linear objective), ``geometric_mean`` and
``arithmetic_mean`` are also available.

For models with a quadratic objective ((MI)QP), the equilibration procedure
differs from the LP case. Instead of scaling the constraint matrix alone, a
symmetric KKT matrix is formed from the quadratic objective matrix and the
constraint matrix, and a modified Ruiz equilibration :cite:`Stellato2020` is
applied jointly to scale the objective and constraints. An additional cost
scaling step normalises the objective coefficients at each iteration.

**Geometric mean** scaling :cite:`Elble2011` scales each row :math:`i` by
:math:`1/\sqrt{a_i^{\max} \cdot a_i^{\min}}`, where :math:`a_i^{\max}` and
:math:`a_i^{\min}` are the largest and smallest absolute nonzero values in that
row or column. This tends to perform well when the coefficient range spans
many orders of magnitude.

**Arithmetic mean** scaling :cite:`Elble2011` scales both rows and columns by
the reciprocal of the mean absolute value of their nonzero entries. It can
perform better than equilibration when the distribution of coefficient
magnitudes is unimodal and roughly symmetric on a log scale.

When a ``geometric_mean`` or ``arithmetic_mean`` method is called for a model
with a quadratic objective, a warning is issued and ``equilibration`` is used
automatically.


Controlling Convergence
***********************

The iterative scaling algorithm repeats until the maximum deviation of the
scaling factors from 1 falls below ``scale_conv_tol``, the number of passes reaches
``scale_passes``, or the elapsed time budget is exhausted. When
``scaling_time_limit`` is set, the algorithm completes the current pass before
stopping; it does not interrupt a pass mid-way.

.. code-block:: python

   m_scaled = gma.scale_model(
       m,
       method="equilibration",
       scale_passes=10,          # allow up to 10 iterations
       scale_conv_tol=1e-6,       # tighter convergence threshold
       scaling_time_limit=30.0,  # stop after completing the current pass
   )

The default of one pass is intentionally conservative: a single scaling pass
can itself take significant time on large models, and the primary goal is to
reduce the solver's work, not to achieve a perfectly scaled matrix at any cost.
The tradeoff has two dimensions:

* **Scaling time vs. solution time.** Additional passes improve conditioning
  and can substantially reduce solver runtime, but they add upfront cost.
  For a model that is already moderately well-scaled, one pass may be enough
  to achieve the desired speedup. For severely ill-conditioned models, allowing
  more passes is likely worthwhile.

* **Conditioning vs. unscaled solution quality.** More scaling passes generally
  improve numerical conditioning but can increase unscaled constraint and bound
  violations after solving, because the solution is recovered in the original
  variable space by applying the inverse scaling transformation.

Users are encouraged to experiment with ``scale_passes`` for their specific
model and to evaluate the resulting ``MaxUnscVio`` alongside solver runtime to
find the right balance.

The ``scaling_lb`` and ``scaling_ub`` parameters bound the scaling factors,
preventing extreme rescaling that could itself introduce numerical issues:

.. code-block:: python

   m_scaled = gma.scale_model(
       m,
       method="geometric_mean",
       scaling_lb=1e-4,   # factors clipped to [1e-4, 1e4]
       scaling_ub=1e4,
   )


Opt-Out Per Variable or Constraint
***********************************

Individual variables or constraints can be excluded from scaling by setting a
``_scale`` attribute to ``0`` before calling :py:func:`~gurobi_modelanalyzer.scale_model`. This is
useful when a particular variable or constraint is already well-scaled and
should not be modified, or when the physical interpretation of a coefficient
must be preserved.

.. code-block:: python

   m.getVarByName("x1")._scale = 0       # exclude from column scaling
   m.getConstrByName("c1")._scale = 0    # exclude from row scaling

   m_scaled = gma.scale_model(m, method="equilibration")

Integer and binary variables are always excluded from column scaling,
regardless of the ``_scale`` attribute. Setting ``_scale = 0`` on a continuous
variable or a constraint simply extends this opt-out to those objects.


User-Provided Initial Scaling
******************************

For cases where domain knowledge about suitable scaling factors is available,
the ``init_scaling`` parameter controls how user-provided initial factors are
used. These factors are set via the ``_init_scaling`` attribute on individual
variables and constraints before calling :py:func:`~gurobi_modelanalyzer.scale_model`.

**Mode 0 (default):** ignore any ``_init_scaling`` attributes and run the
iterative algorithm from the identity scaling.

.. code-block:: python

   m_scaled = gma.scale_model(m, method="equilibration", init_scaling=0)

**Mode 1:** apply ``_init_scaling`` as the final scaling and return
immediately, without running any iterative algorithm. This is useful when
the user has pre-computed scaling factors and wants to apply them directly.

.. code-block:: python

   for var in m.getVars():
       var._init_scaling = my_col_factors[var.VarName]
   for constr in m.getConstrs():
       constr._init_scaling = my_row_factors[constr.ConstrName]

   m_scaled = gma.scale_model(m, method="equilibration", init_scaling=1)

**Mode 2 (warmstart):** pre-apply ``_init_scaling`` to the coefficient matrix,
then run the iterative algorithm on top. The final scaling factors are the
product of the user-provided values and the algorithm's output. This is useful
when rough scaling estimates are available but further refinement is desired.

.. code-block:: python

   m_scaled = gma.scale_model(m, method="equilibration", init_scaling=2)

When both ``_init_scaling`` and ``_scale = 0`` are set on the same variable or
constraint, ``_init_scaling`` takes priority. The scaling factor is fixed at
the value of ``_init_scaling`` and is not modified by the algorithm. It is
held constant throughout all passes, effectively locking that factor in place.


Accessing the Scaling Matrices
*******************************

The row and column scaling factors are stored as diagonal sparse matrices on
the :ref:`ScaledModel <APIScaledModelLabel>` object and can be accessed via
the ``ColScaling`` and ``RowScaling`` properties:

.. code-block:: python

   m_scaled = gma.scale_model(m, method="geometric_mean")

   col_factors = m_scaled.ColScaling.diagonal()   # shape: (num_vars,)
   row_factors = m_scaled.RowScaling.diagonal()   # shape: (num_constrs,)

   print(f"Column factor range: [{col_factors.min():.2e}, {col_factors.max():.2e}]")
   print(f"Row factor range:    [{row_factors.min():.2e}, {row_factors.max():.2e}]")

A wide column factor range suggests the original model had highly variable
coefficient magnitudes across variables. Inspecting the individual factors can
reveal which variables or constraints drove the need for scaling.

The factor for each individual variable or constraint is also accessible
directly as a ``scaling_factor`` attribute on the wrapper objects returned by
:py:meth:`ScaledModel.getVarsUnscaled` and
:py:meth:`ScaledModel.getConstrsUnscaled`:

.. code-block:: python

   m_scaled = gma.scale_model(m, method="equilibration")

   for var in m_scaled.getVarsUnscaled():
       print(f"{var.VarName}: scaling factor = {var.scaling_factor:.4e}")

   for constr in m_scaled.getConstrsUnscaled():
       print(f"{constr.ConstrName}: scaling factor = {constr.scaling_factor:.4e}")

   # Quadratic constraints follow the same pattern:
   for qconstr in m_scaled.getQConstrsUnscaled():
       print(f"{qconstr.QCName}: scaling factor = {qconstr.scaling_factor:.4e}")


Computing the Unscaled Objective
*********************************

To retrieve the objective value in the original variable space after
optimization, use :py:meth:`ScaledModel.computeUnscObj`:

.. code-block:: python

   m_scaled.optimize()
   m_scaled.computeUnscObj()
   print(f"Unscaled objective: {m_scaled.UnscObjVal:.6e}")
   print(f"Scaled objective:   {m_scaled.ObjVal:.6e}")

The unscaled and scaled objective values will generally differ because the
objective coefficients are also transformed during scaling.


Scaling Log
***********

Scaling progress is printed to the console by default
(``scaling_log_to_console=1``). It can be redirected to a file, or suppressed
entirely:

.. code-block:: python

   # Log to file only
   m_scaled = gma.scale_model(
       m,
       method="equilibration",
       scaling_log="scaling.log",
       scaling_log_to_console=0,
   )

   # Log to both console and file
   m_scaled = gma.scale_model(
       m,
       method="equilibration",
       scaling_log="scaling.log",
       scaling_log_to_console=1,
   )

The log reports the original model's coefficient ranges, the maximum deviation
of the scaling factors from 1 per pass, elapsed time per pass, total scaling
time, and the scaled model's coefficient ranges. A typical log looks like::

   Scaling Method: arithmetic_mean
   Scale Passes:   5
   Conv. Tol.:     1.000000e-04
   Original Model Statistics:
   Statistics for model 'glass4':
     Problem type                : MIP
     Linear constraint matrix    : 396 rows, 322 columns, 1815 nonzeros
     Variable types              : 20 continuous, 302 integer (0 binary)
     Matrix range                : [1e+00, 8e+06]
     Objective range             : [1e+00, 1e+06]
     Bounds range                : [1e+00, 8e+02]
     RHS range                   : [1e+00, 8e+06]
   Scale Pass   Max Factor Dev.   Time (s)
   1            8.104301e+03      0.01
   2            3.270924e+01      0.01
   3            2.626976e+01      0.02
   4            4.890177e+00      0.02
   5            1.053361e+00      0.03
   Building scaled model...
   Scaling completed in 0.04 seconds
   Scaled Model Ranges:
     Matrix range                : [2e-01, 4e+00]
     Objective range             : [2e+04, 1e+10]
     Bounds range                : [5e-03, 1e+00]
     RHS range                   : [9e-03, 4e+00]

The **Max Factor Dev.** column shows the maximum deviation of the scaling
factors from 1 in that pass. When it falls below ``scale_conv_tol``,
the algorithm has converged and no further passes are performed. The **Time (s)**
column shows cumulative wall-clock time up to and including that pass. Comparing the
original and scaled **Matrix range** shows how much the coefficient spread has
been reduced. A tighter range often indicates better conditioning.


.. _ScalingFilesLabel:

Scaling Files
*************

Scaling factors can be saved to and loaded from plain-text ``.scl`` files.
This makes it possible to:

* reproduce a previously found scaling exactly on a re-run,
* share scaling factors between runs or users,
* provide domain-knowledge-based initial factors without writing Python code,
* use scaling factors found via the Python API as input to ``gurobi_cls``,
  and vice versa.


File Format
~~~~~~~~~~~

A ``.scl`` file is a plain-text file. Lines beginning with ``#`` and blank
lines are ignored. An optional version header must appear before the first
section if present::

   GRB_SCL_FILE_VERSION 1

Data is organised into up to three named sections::

   SECTION VARS
   SECTION CONSTRS
   SECTION QCONSTRS

Each data line within a section has the form::

   name  factor  lock_flag

where ``name`` is the variable or constraint name as it appears in the model,
``factor`` is a positive floating-point scaling factor, and ``lock_flag`` is
either ``0`` (keep the factor fixed; the algorithm will not modify it) or
``1`` (use the factor as a warmstart; the algorithm may adjust it further).

All three sections are optional. A file may contain only a subset of the
variables or constraints in the model; any object not listed defaults to an
initial factor of 1.0.

A minimal example::

   # My custom scaling factors
   GRB_SCL_FILE_VERSION 1

   SECTION VARS
   price    1e-3  0
   quantity 1e+2  1

   SECTION CONSTRS
   budget   5e-1  0


Writing Scaling Factors (Python API)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After calling :py:func:`~gurobi_modelanalyzer.scale_model`, use
:py:meth:`ScaledModel.write_scaling` to export the computed factors:

.. code-block:: python

   import gurobipy as gp
   import gurobi_modelanalyzer as gma

   m = gp.read("model.mps")
   m_scaled = gma.scale_model(m, method="equilibration")

   # Save with lock_flag=0 (default): factors are fixed on re-import
   m_scaled.write_scaling("model.scl")

   # Save with lock_flag=1: factors act as a warmstart on re-import
   m_scaled.write_scaling("model.scl", lock_factors=False)


Reading Scaling Factors (Python API)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use :py:func:`~gurobi_modelanalyzer.read_scaling_file` to load a ``.scl``
file and apply its factors to a model before calling
:py:func:`~gurobi_modelanalyzer.scale_model`:

.. code-block:: python

   import gurobipy as gp
   from gurobi_modelanalyzer.scaling import scale_model, read_scaling_file

   m = gp.read("model.mps")
   read_scaling_file("model.scl", m)
   m_scaled = scale_model(m, method="equilibration", init_scaling=2)

The function sets ``_init_scaling`` and, where applicable, ``_scale = 0``
directly on the model's variable and constraint objects. The ``init_scaling``
parameter of :py:func:`~gurobi_modelanalyzer.scale_model` must be set to ``1``
or ``2`` for these attributes to take effect (see
`User-Provided Initial Scaling`_ above). Use ``init_scaling=2`` to run the
algorithm as a warmstart on top of the loaded factors, or ``init_scaling=1``
to apply them without any further iteration.

Malformed lines and unrecognised names issue :class:`UserWarning` automatically
via Python's ``warnings`` module. They can be suppressed with::

   import warnings
   warnings.filterwarnings("ignore", category=UserWarning)


Round-Trip Example
~~~~~~~~~~~~~~~~~~

The following pattern scales a model, saves the result, and reproduces the
same scaling on a later run:

.. code-block:: python

   import gurobipy as gp
   import gurobi_modelanalyzer as gma
   from gurobi_modelanalyzer.scaling import read_scaling_file

   # --- First run: compute and save scaling ---
   m = gp.read("model.mps")
   m_scaled = gma.scale_model(m, method="equilibration")
   m_scaled.write_scaling("model.scl")          # lock_flag=0 by default

   # --- Later run: reproduce the exact same scaling ---
   m2 = gp.read("model.mps")
   read_scaling_file("model.scl", m2)
   m2_scaled = gma.scale_model(m2, method="equilibration", init_scaling=2)

Because the saved factors use ``lock_flag=0``, the algorithm cannot modify
them further, and ``m2_scaled`` has the same coefficient matrix as ``m_scaled``.
