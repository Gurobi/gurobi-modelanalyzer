.. _ScalingQSGuideLabel:

Quick Start Guide
#################

.. _ScalingQSIntroductionLabel:

Introduction
************

The Model Scaling module can be used to scale the coefficient matrix of a
Gurobi model to improve numerical conditioning prior to solving. Poorly scaled
models — those with large differences in coefficient magnitudes across rows or
columns — can cause numerical instability in the solver, leading to inaccurate
solutions or convergence issues.

The module works by computing diagonal row and column scaling matrices
:math:`D_r` and :math:`D_c` such that the transformed constraint matrix
:math:`\hat{A} = D_r A D_c` has improved coefficient ranges. The resulting
scaled model is returned as a :ref:`ScaledModel <APIScaledModelLabel>` object,
which is a subclass of ``gurobipy.Model`` and can be optimized directly. After
solving, the original (unscaled) solution can be recovered via
:py:meth:`ScaledModel.getVarsUnscaled`.

Three scaling methods are available:

* **equilibration**: iteratively scales rows and columns to bring coefficient
  magnitudes to a similar range. Supports (MI)LP, (MI)QP, and (MI)QCP models.
  This is the recommended starting point and the only method supported for
  models with a quadratic objective.
* **geometric_mean**: scales rows and columns by the reciprocal of the
  geometric mean of their minimum and maximum absolute values.
  Supports (MI)LP and (MI)QCP models.
* **arithmetic_mean**: scales both rows and columns by the reciprocal of
  their mean absolute value. Supports (MI)LP and (MI)QCP models.

See the :doc:`advanced_usage_scaling` section for guidance on choosing between
methods and controlling the scaling algorithm in detail.


Scaling a Model
***************

The simplest way to use the module is to call :py:func:`~gurobi_modelanalyzer.scale_model` with a
model and a method name:

.. code-block:: python

   import gurobipy as gp
   import gurobi_modelanalyzer as gma

   m = gp.read("mymodel.mps")
   m_scaled = gma.scale_model(m, method="equilibration")

   m_scaled.optimize()

:py:func:`~gurobi_modelanalyzer.scale_model` returns a fully constructed
:ref:`ScaledModel <APIScaledModelLabel>` object. You can inspect it before
solving — for example, examine the scaling factors, check coefficient ranges,
or inspect variables and constraints — and call ``optimize()`` only when ready.

The scaled model is solved in the transformed space. To retrieve the solution
values in the original variable space:

.. code-block:: python

   for var in m_scaled.getVarsUnscaled():
       print(f"{var.VarName}: unscaled = {var.Xunsc:.6e}, scaled = {var.X:.6e}")


Checking Solution Quality
*************************

After optimization, use :py:meth:`ScaledModel.computeUnscVio` to compute
constraint and bound violations in the original (unscaled) variable space:

.. code-block:: python

   m_scaled.computeUnscVio()

   print(f"Max constraint violation: {m_scaled.MaxUnscConstrVio:.2e}")
   print(f"Max bound violation:      {m_scaled.MaxUnscBoundVio:.2e}")

``MaxUnscVio`` is simply the larger of the two — equivalent to
``max(MaxUnscConstrVio, MaxUnscBoundVio)``.

Individual violations per constraint are accessible via the wrappers returned
by :py:meth:`ScaledModel.getConstrsUnscaled`:

.. code-block:: python

   for c in m_scaled.getConstrsUnscaled():
       if c.UnscViolation > 1e-6:
           print(f"{c.ConstrName}: violation = {c.UnscViolation:.2e}")
