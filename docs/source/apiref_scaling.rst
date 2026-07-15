.. _ScalingAPIRefLabel:

API Reference
#############


.. _APIscale_modelLabel:

.. py:function:: gurobi_modelanalyzer.scale_model(model, method, scale_passes=1, scale_conv_tol=1e-4, scaling_lb=1e-8, scaling_ub=1e8, value_threshold=1e-13, scaling_time_limit=inf, scaling_log="", scaling_log_to_console=1, init_scaling=0, env=None)

   Scale a Gurobi optimization model to improve numerical conditioning.

   Creates a scaled copy of the input model using the specified scaling
   method. The scaled model can be optimized directly and provides methods
   to recover the solution in the original (unscaled) variable space.

   :param model: Required Gurobi model to scale.
   :param method: Scaling method to use. One of:

                  ``'equilibration'``: iteratively scales rows and columns to bring
                  coefficient magnitudes to a similar range. Supports all model
                  types: (MI)LP, (MI)QP, (MI)QCP, and (MI)QCQP.

                  ``'geometric_mean'``: scales rows and columns by
                  :math:`1/\sqrt{a_i^{\max} \cdot a_i^{\min}}`, where
                  :math:`a_i^{\max}` and :math:`a_i^{\min}` are the largest
                  and smallest absolute nonzero values in the row or column.
                  Supports (MI)LP and (MI)QCP models (linear objective only).

                  ``'arithmetic_mean'``: scales rows and columns by the
                  reciprocal of the mean absolute value of their nonzero entries.
                  Supports (MI)LP and (MI)QCP models (linear objective only).

                  See :ref:`ScalingAdvUsageLabel` for the model class taxonomy
                  and a comparison of methods.

   :param scale_passes: Maximum number of scaling iterations. Default: 1.
   :param scale_conv_tol: Convergence tolerance. Scaling stops early when the
                         maximum deviation of the scaling factors from 1 falls
                         below this threshold. Default: 1e-4.
   :param scaling_lb: Lower bound for scaling factors. Prevents extreme
                      downscaling. Default: 1e-8.
   :param scaling_ub: Upper bound for scaling factors. Prevents extreme
                      upscaling. Default: 1e8.
   :param value_threshold: Coefficients with absolute value below this
                           threshold are treated as zero. Default: 1e-13.
   :param scaling_time_limit: Time limit in seconds for the scaling
                              iterations. If reached, the best scaling found
                              so far is used. Default: no limit.
   :param scaling_log: Optional path to a log file. If provided, scaling
                       progress is written to this file. Default: no file.
   :param scaling_log_to_console: Set to 1 (default) to print scaling
                                  progress to the console, 0 to suppress.
   :param init_scaling: Controls use of user-provided initial scaling factors
                        set via the ``_init_scaling`` attribute on variables
                        and constraints:

                        ``0`` (default): ignore ``_init_scaling``; run the
                        iterative algorithm from the identity scaling.

                        ``1``: apply ``_init_scaling`` as the final scaling
                        and return immediately without running the iterative
                        algorithm.

                        ``2`` (warmstart): pre-apply ``_init_scaling``, then
                        run the iterative algorithm on top. The final factors
                        are the product of the user-provided values and the
                        algorithm's output.

   :param env: Optional Gurobi environment (``gurobipy.Env``) to use for
               the scaled model.
   :return: A :ref:`ScaledModel <APIScaledModelLabel>` object containing the
            scaled model with scaling information attached.


.. _APIread_scaling_fileLabel:

.. py:function:: gurobi_modelanalyzer.read_scaling_file(path, model)

   Parse a ``.scl`` scaling input file and apply the specified initial scaling
   factors to the Gurobi model's variables and constraints.

   The function sets ``_init_scaling`` and, where applicable, ``_scale``
   attributes directly on the ``gurobipy.Var`` and ``gurobipy.Constr`` objects
   of *model*. After calling this function, pass ``init_scaling=2`` to
   :py:func:`~gurobi_modelanalyzer.scale_model` to run the iterative algorithm
   as a warmstart on top of the loaded factors, or ``init_scaling=1`` to apply
   them without any further iteration.

   Malformed lines and unrecognised variable or constraint names issue
   :class:`UserWarning` via Python's standard ``warnings`` machinery.

   See :ref:`ScalingFilesLabel` for a description of the ``.scl`` file format
   and usage examples.

   :param path: Path to the ``.scl`` scaling input file.
   :type path: str
   :param model: Gurobi model whose objects will receive ``_init_scaling`` /
                 ``_scale`` attributes.
   :return: ``None``


.. _APIScaledModelLabel:

ScaledModel
***********

``ScaledModel`` is a subclass of ``gurobipy.Model`` returned by
:py:func:`~gurobi_modelanalyzer.scale_model`. It adds methods and properties for recovering
unscaled solutions and computing violations in the original variable space.

.. py:method:: ScaledModel.getVarsUnscaled()

   Return a list of :ref:`ScaledVar <APIScaledVarLabel>` objects, one per
   variable in the model. Each object exposes both the scaled solution value
   (``X``) and the unscaled value (``Xunsc``). Must be called after
   optimization.

   :return: List of :ref:`ScaledVar <APIScaledVarLabel>` objects.

.. py:method:: ScaledModel.getConstrsUnscaled()

   Return a list of :ref:`ScaledConstr <APIScaledConstrLabel>` objects, one
   per linear constraint. After calling :py:meth:`computeUnscVio`, each
   object exposes the unscaled constraint violation via ``UnscViolation``.

   :return: List of :ref:`ScaledConstr <APIScaledConstrLabel>` objects.

.. py:method:: ScaledModel.getQConstrsUnscaled()

   Return a list of :ref:`ScaledQConstr <APIScaledConstrLabel>` objects, one
   per quadratic constraint. After calling :py:meth:`computeUnscVio`, each
   object exposes the unscaled constraint violation via ``UnscViolation``.

   :return: List of :ref:`ScaledQConstr <APIScaledConstrLabel>` objects.

.. py:method:: ScaledModel.computeUnscVio()

   Compute constraint and bound violations in the original (unscaled) variable
   space. Must be called after optimization. Populates ``UnscViolation`` on
   all constraint wrappers and ``UnscBoundViolation`` on all variable wrappers,
   and sets the ``MaxUnscVio``, ``MaxUnscConstrVio``, and ``MaxUnscBoundVio``
   properties. The original model is stored automatically by
   :func:`~gurobi_modelanalyzer.scale_model`.

.. py:method:: ScaledModel.computeUnscObj()

   Compute the objective value in the original (unscaled) variable space using
   the unscaled solution values from :py:meth:`getVarsUnscaled`. Must be
   called after optimization. Result is stored in :py:attr:`ScaledModel.UnscObjVal`.

   :return: ``None`` (access the result via :py:attr:`ScaledModel.UnscObjVal`).

.. py:method:: ScaledModel.write_scaling(path, lock_factors=True)

   Export the scaling factors computed by :py:func:`~gurobi_modelanalyzer.scale_model`
   to a ``.scl`` file. The file can later be passed to
   :py:func:`~gurobi_modelanalyzer.read_scaling_file` or to ``gurobi_cls`` via
   ``--scaling-file`` to reproduce or continue from the same scaling.

   All variables and constraints are written, including those with a factor of
   1.0. When ``lock_factors=True``, a factor of 1.0 with ``lock_flag=0`` locks
   that object at the identity scaling and prevents the algorithm from modifying
   it on re-import.

   :param path: Output file path (conventionally with ``.scl`` extension).
   :type path: str
   :param lock_factors: If ``True`` (default), all entries are written with
                        ``lock_flag=0``, meaning the factors are kept fixed
                        when the file is re-imported and the algorithm cannot
                        modify them. If ``False``, ``lock_flag=1`` is written
                        so the factors act as a warmstart.
   :type lock_factors: bool
   :return: ``None``

.. py:attribute:: ScaledModel.UnscObjVal

   Unscaled objective value computed by :py:meth:`computeUnscObj`. ``None``
   until that method is called.

.. py:attribute:: ScaledModel.MaxUnscVio

   Maximum unscaled violation across all constraints and variable bounds.
   Available after calling :py:meth:`computeUnscVio`.

.. py:attribute:: ScaledModel.MaxUnscConstrVio

   Maximum unscaled violation across all linear and quadratic constraints.
   Available after calling :py:meth:`computeUnscVio`.

.. py:attribute:: ScaledModel.MaxUnscBoundVio

   Maximum unscaled variable bound violation.
   Available after calling :py:meth:`computeUnscVio`.

.. py:attribute:: ScaledModel.ScalingTime

   Wall-clock time in seconds taken by the scaling procedure.

.. py:attribute:: ScaledModel.ColScaling

   Diagonal column scaling matrix as a ``scipy.sparse`` matrix. Entry
   :math:`i` contains the scaling factor applied to variable :math:`i`.

.. py:attribute:: ScaledModel.RowScaling

   Diagonal row scaling matrix as a ``scipy.sparse`` matrix. Entry
   :math:`i` contains the scaling factor applied to constraint :math:`i`.


.. _APIScaledVarLabel:

ScaledVar
*********

Wrapper around a ``gurobipy.Var`` object returned by
:py:meth:`ScaledModel.getVarsUnscaled`. All standard Gurobi variable
attributes (e.g. ``VarName``, ``LB``, ``UB``) are forwarded to the
underlying variable.

.. py:attribute:: ScaledVar.X

   Solution value in the scaled model space.

.. py:attribute:: ScaledVar.Xunsc

   Solution value recovered in the original (unscaled) space:
   :math:`x_i = s_i \cdot y_i`, where :math:`s_i` is the column scaling
   factor and :math:`y_i` is the scaled solution value.

.. py:attribute:: ScaledVar.scaling_factor

   Column scaling factor :math:`s_i` applied to this variable by
   :func:`~gurobi_modelanalyzer.scale_model`. Always positive.

.. py:attribute:: ScaledVar.UnscBoundViolation

   Unscaled bound violation for this variable. Available after calling
   :py:meth:`ScaledModel.computeUnscVio`.


.. _APIScaledConstrLabel:

ScaledConstr / ScaledQConstr
****************************

Wrappers around ``gurobipy.Constr`` and ``gurobipy.QConstr`` objects returned
by :py:meth:`ScaledModel.getConstrsUnscaled` and
:py:meth:`ScaledModel.getQConstrsUnscaled` respectively. All standard
Gurobi constraint attributes are forwarded to the underlying object.

.. py:attribute:: ScaledConstr.UnscViolation
                  ScaledQConstr.UnscViolation

   Unscaled constraint violation. Available after calling
   :py:meth:`ScaledModel.computeUnscVio`.

.. py:attribute:: ScaledConstr.scaling_factor
                  ScaledQConstr.scaling_factor

   Row scaling factor applied to this constraint by
   :func:`~gurobi_modelanalyzer.scale_model`. Always positive.
