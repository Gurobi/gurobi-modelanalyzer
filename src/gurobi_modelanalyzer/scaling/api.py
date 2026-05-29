import gurobipy as gp
from gurobipy import GRB
import scipy.sparse
import numpy as np
from typing import List, Tuple, Dict, Union
import warnings
import time
import io
import logging
import sys

import gurobi_modelanalyzer.common as common

from .methods import (
    ModelData,
    _threshold_small_coefficients,
    _iterative_scaling,
    quad_equilibration,
    _print_scaling_log,
    _extract_range_stats,
    _scale_single_qconstr,
)
from .scaled_wrappers import ScaledModel

logger = logging.getLogger("gurobi_modelanalyzer.scaling")
logger.setLevel(logging.DEBUG)
logger.propagate = False


def _build_init_scaling_matrices(
    gurobi_vars: List[gp.Var],
    gurobi_constrs: List[gp.Constr],
) -> Tuple[scipy.sparse.dia_matrix, scipy.sparse.dia_matrix]:
    """
    Build diagonal initial scaling matrices from ``_init_scaling``
    attributes set on Gurobi variable and constraint objects.

    Entries without ``_init_scaling`` default to 1.0 (no effect).
    """
    col_factors = np.array(
        [float(getattr(v, "_init_scaling", 1.0)) for v in gurobi_vars]
    )
    row_factors = np.array(
        [float(getattr(c, "_init_scaling", 1.0)) for c in gurobi_constrs]
    )
    return scipy.sparse.diags(col_factors), scipy.sparse.diags(row_factors)


def _qconstr_row_scale_override(
    qconstr: gp.QConstr,
    init_scaling: int,
) -> float:
    """
    Return the row_scale_override value for a quadratic constraint, or
    None if the normal (norm-based) computation should be used.

    Rules:
    - init_scaling == 1: always use _init_scaling when set
      (user-provided scaling only, no algorithm).
    - init_scaling == 2 with _scale == 0: use _init_scaling to lock the
      row factor (the variable is excluded from further row scaling).
    - All other cases: return None (norm-based or skip_row_scale applies).
    """
    init_val = getattr(qconstr, "_init_scaling", None)
    if init_val is None:
        return None
    if init_scaling == 1:
        return float(init_val)
    if init_scaling == 2 and getattr(qconstr, "_scale", 1) == 0:
        return float(init_val)
    return None


def _capture_model_stats(model: gp.Model) -> str:
    """
    Capture the output of model.printStats() as a string.

    Parameters:
    -----------
    model : gp.Model
        The Gurobi model to analyze

    Returns:
    --------
    str
        The statistics output as a string
    """
    # Ensure model is updated before capturing stats
    model.update()

    # Capture stdout
    old_stdout = sys.stdout
    sys.stdout = buffer = io.StringIO()

    try:
        model.printStats()
        output = buffer.getvalue()
    finally:
        sys.stdout = old_stdout

    return output


def scale_model(
    model: gp.Model,
    method: str,
    scale_passes: int = 1,
    scale_conv_tol: float = 1e-4,
    scaling_lb: float = 1e-8,
    scaling_ub: float = 1e8,
    value_threshold: float = 1e-13,
    scaling_time_limit: float = float("inf"),
    scaling_log: str = "",
    scaling_log_to_console: int = 1,
    init_scaling: int = 0,
    env: gp.Env = None,
) -> ScaledModel:
    """
    Scale a Gurobi optimization model to improve numerical conditioning.

    Creates a scaled version of the input model using the specified
    scaling method. The scaled model can be solved, and the solution
    can be unscaled back to the original variable space.

    Parameters:
    -----------
    model : gp.Model
        The Gurobi model to scale
    method : str
        Scaling method to use. Options:
        - 'equilibration': Mean-based equilibration (works for LP, QP, QCP)
        - 'geometric_mean': Geometric mean scaling (LP, QCP; not QP)
        - 'arithmetic_mean': Arithmetic mean scaling (LP, QCP; not QP)
    scale_passes : int, optional
        Maximum number of scaling iterations (default: 5)
    scale_conv_tol : float, optional
        Convergence tolerance: scaling stops early when the maximum deviation
        of the scaling factors from 1 falls below this threshold (default: 1e-4)
    scaling_lb : float, optional
        Lower bound for scaling factors to avoid extreme values (default: 1e-8)
    scaling_ub : float, optional
        Upper bound for scaling factors to avoid extreme values (default: 1e8)
    value_threshold : float, optional
        Threshold below which coefficients are set to zero (default: 1e-13)
    scaling_time_limit : float, optional
        Time limit in seconds for scaling iterations. Scaling will
        stop when this limit is reached and use the latest scaled
        matrices (default: inf - no limit)
    scaling_log : str, optional
        File path to write scaling log to. Empty string means no
        file output (default: "")
    scaling_log_to_console : int, optional
        1 to print scaling log to console, 0 to suppress console
        output (default: 1)
    init_scaling : int, optional
        Controls how user-provided initial scaling factors are used:
        - 0 (default): ignore any _init_scaling attributes; run the
          algorithm from the identity scaling.
        - 1: read _init_scaling from variables and constraints, apply
          that scaling, and return immediately without running any
          iterative algorithm.
        - 2 (warmstart): read _init_scaling, pre-apply that scaling to
          the matrices, then run the iterative algorithm on top.
          The final scaling stored in the model is the product of the
          user-provided initial scaling and the algorithm scaling.
    env : gp.Env, optional
        Gurobi environment to use for the scaled model.

    Returns:
    --------
    ScaledModel
        A scaled version of the input model with scaling information attached.
        Use getVarsUnscaled() to access unscaled solution values.

    Notes:
    ------
    - For models with quadratic objectives, only 'equilibration'
      method is supported
    - For models with quadratic constraints (but no quadratic
      objective), all methods work
    - Integer and binary variables are never scaled regardless of any
      _scale attribute
    - Set var._scale = 0 on a continuous variable to exclude it from
      column scaling
    - Set constr._scale = 0 on a linear or quadratic constraint to
      exclude it from row scaling
    - Set var._init_scaling = k on a variable to provide an initial
      column scaling factor of k (used when init_scaling > 0)
    - Set constr._init_scaling = k on a linear or quadratic constraint
      to provide an initial row scaling factor of k (used when
      init_scaling > 0)
    - _init_scaling takes priority over _scale=0: a variable/constraint
      with both set will use _init_scaling as its (locked) factor
    - For QCP constraints in mode 2, _init_scaling only locks the row
      factor when _scale=0 is also set; otherwise the norm-based
      computation runs on the column-pre-scaled constraint
    - The scaled model includes scaling matrices stored as
      _col_scaling and _row_scaling attributes
    """
    if init_scaling not in (0, 1, 2):
        raise ValueError(f"init_scaling must be 0, 1, or 2 (got {init_scaling!r}).")

    # Start timing
    total_start_time = time.time()

    # Ensure model is updated
    model.update()

    # Manage OutputFlag: temporarily enable if needed for logging
    original_output_flag = model.Params.OutputFlag
    needs_output = scaling_log_to_console or scaling_log
    if original_output_flag == 0 and needs_output:
        model.setParam("OutputFlag", 1)

    # Configure logging handlers for this scaling call
    _log_handlers: List[logging.Handler] = []
    if scaling_log_to_console:
        _console_handler = logging.StreamHandler()
        _console_handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(_console_handler)
        _log_handlers.append(_console_handler)
    if scaling_log:
        _file_handler = logging.FileHandler(scaling_log, mode="w")
        _file_handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(_file_handler)
        _log_handlers.append(_file_handler)

    # Get model data from original model
    model_data = ModelData.from_gurobi_model(model)

    # Capture original model statistics
    original_stats = _capture_model_stats(model)

    # Compute columns to scale: skip integer/binary unconditionally;
    # skip continuous variables explicitly marked with _scale=0.
    gurobi_vars = model.getVars()
    cols_to_scale = []
    for i, var_type in enumerate(model_data.var_types):
        if var_type in (GRB.INTEGER, GRB.BINARY):
            continue
        if getattr(gurobi_vars[i], "_scale", 1) == 0:
            continue
        cols_to_scale.append(i)

    # Compute rows to scale: skip any constraint marked with _scale=0.
    rows_to_scale = []
    for i, constr in enumerate(model.getConstrs()):
        if getattr(constr, "_scale", 1) == 0:
            continue
        rows_to_scale.append(i)

    # Build initial scaling matrices from _init_scaling attributes when needed.
    if init_scaling > 0:
        col_init_scaling, row_init_scaling = _build_init_scaling_matrices(
            gurobi_vars, model.getConstrs()
        )
    else:
        col_init_scaling = row_init_scaling = None

    # Check if model has quadratic objective terms
    q_matrix = model.getQ()

    iteration_logs = []

    # Print log header if logging is enabled
    _print_scaling_log(
        method,
        original_stats,
        scale_passes,
        [],
        0.0,
        "",
        scaling_time_limit,
        scale_conv_tol,
        mode="header",
    )

    if init_scaling == 1:
        # Mode 1: apply user-provided scaling only, skip the algorithm.
        logger.info("Using user-provided initial scaling (init_scaling=1).")
        col_scaling = col_init_scaling
        row_scaling = row_init_scaling
        scaled_matrix = row_scaling @ model_data.constr_matrix @ col_scaling
        obj_vector_scaled = col_scaling @ model_data.obj_vector
        if q_matrix.nnz > 0:
            scaled_q_matrix = col_scaling @ q_matrix @ col_scaling
    else:
        # Modes 0 and 2: run the iterative algorithm.
        # Mode 2 (warmstart): pre-scale input matrices with user factors so
        # the algorithm continues from that starting point.
        if init_scaling == 2:
            constr_matrix_in = (
                row_init_scaling @ model_data.constr_matrix @ col_init_scaling
            )
            if q_matrix.nnz > 0:
                q_matrix_in = col_init_scaling @ q_matrix @ col_init_scaling
                obj_vector_in = col_init_scaling @ model_data.obj_vector
            else:
                q_matrix_in = q_matrix
                obj_vector_in = model_data.obj_vector
        else:  # mode 0: use original matrices unchanged
            constr_matrix_in = model_data.constr_matrix
            q_matrix_in = q_matrix
            obj_vector_in = model_data.obj_vector

        if q_matrix.nnz == 0:  # LP / QCP: no quadratic objective
            (scaled_matrix, row_scaling, col_scaling, iteration_logs) = (
                _iterative_scaling(
                    constr_matrix_in,
                    cols_to_scale,
                    rows_to_scale,
                    scale_passes,
                    scale_conv_tol,
                    method,
                    scaling_time_limit=scaling_time_limit,
                )
            )
            # Accumulate initial factors into the total for mode 2
            if init_scaling == 2:
                col_scaling = col_init_scaling @ col_scaling
                row_scaling = row_scaling @ row_init_scaling
            obj_vector_scaled = col_scaling @ model_data.obj_vector
        else:  # QP: quadratic objective
            if method != "equilibration":
                warnings.warn(
                    "Equilibration is the only supported method for "
                    "quadratic objectives. Using equilibration instead.",
                    UserWarning,
                )
            (
                scaled_matrix,
                scaled_q_matrix,
                obj_vector_scaled,
                row_scaling,
                col_scaling,
                iteration_logs,
            ) = quad_equilibration(
                constr_matrix_in,
                obj_vector_in,
                q_matrix_in,
                cols_to_scale,
                rows_to_scale,
                scale_passes,
                scale_conv_tol,
                scaling_lb=scaling_lb,
                scaling_ub=scaling_ub,
                scaling_time_limit=scaling_time_limit,
            )
            # Accumulate initial factors into the total for mode 2
            if init_scaling == 2:
                col_scaling = col_init_scaling @ col_scaling
                row_scaling = row_scaling @ row_init_scaling

    # Print separator before model building phase
    logger.info("Building scaled model...")

    # Compute scaled data
    rhs_vector_scaled = row_scaling @ model_data.rhs_vector
    col_diag = col_scaling.diagonal()
    # Avoid extreme inverse values
    col_diag_safe = np.clip(col_diag, scaling_lb, scaling_ub)
    col_scaling_inv = scipy.sparse.diags(1.0 / col_diag_safe)

    lb_vector_scaled = col_scaling_inv @ model_data.lb_vector
    ub_vector_scaled = col_scaling_inv @ model_data.ub_vector
    var_names_scaled = [name + "_scaled" for name in model_data.var_names]
    constr_names_scaled = [name + "_scaled" for name in model_data.constr_names]

    # Clean small coefficients
    scaled_matrix = _threshold_small_coefficients(scaled_matrix, value_threshold)
    rhs_vector_scaled = _threshold_small_coefficients(
        rhs_vector_scaled, value_threshold
    )
    obj_vector_scaled = _threshold_small_coefficients(
        obj_vector_scaled, value_threshold
    )
    lb_vector_scaled = _threshold_small_coefficients(lb_vector_scaled, value_threshold)
    ub_vector_scaled = _threshold_small_coefficients(ub_vector_scaled, value_threshold)
    if q_matrix.nnz > 0:
        scaled_q_matrix = _threshold_small_coefficients(
            scaled_q_matrix, value_threshold
        )

    # Create linear terms of ScaledModel with scaled data using matrix API
    model_scaled = ScaledModel(model.ModelName + "_scaled", env=env)

    # Add variables with scaled bounds and objective
    vars_list = []
    for i in range(len(model_data.var_names)):
        var = model_scaled.addVar(
            lb=lb_vector_scaled[i],
            ub=ub_vector_scaled[i],
            obj=obj_vector_scaled[i],
            vtype=model_data.var_types[i],
            name=var_names_scaled[i],
        )
        vars_list.append(var)

    # Set objective sense
    model_scaled.ModelSense = model.ModelSense
    model_scaled.update()

    # Add scaled constraints using matrix API
    model_scaled.addMConstr(
        scaled_matrix,
        vars_list,
        model_data.constr_sense,
        rhs_vector_scaled,
        constr_names_scaled,
    )

    model_scaled.update()

    # Store scaling matrices and reference to the original model
    model_scaled._col_scaling = col_scaling
    model_scaled._row_scaling = row_scaling
    model_scaled._original_model = model

    # Add quadratic objective term if present
    if q_matrix.nnz > 0:
        lin_objective = model_scaled.getObjective()
        # Get matrix variable representation
        x_mvars = gp.MVar(model_scaled.getVars())
        full_objective = lin_objective + x_mvars.T @ scaled_q_matrix @ x_mvars
        model_scaled.setObjective(full_objective, model.ModelSense)
        model_scaled.update()
    # Scale quadratic constraints if present
    if model.isQCP:
        logger.info("Scaling quadratic constraints...")
        qconstrs = model.getQConstrs()
        # Process quadratic constraints
        qconstr_results = [
            _scale_single_qconstr(
                qconstr,
                model,
                col_scaling,
                skip_row_scale=(getattr(qconstr, "_scale", 1) == 0),
                row_scale_override=_qconstr_row_scale_override(qconstr, init_scaling),
                scaling_lb=scaling_lb,
            )
            for qconstr in qconstrs
        ]

        # Add scaled quadratic constraints to model
        quad_scaling_factors = []
        for (
            qc_scaled,
            q_scaled,
            sense,
            rhs_scaled,
            scaling_factor,
            name,
        ) in qconstr_results:
            model_scaled.addMQConstr(qc_scaled, q_scaled, sense, rhs_scaled, name=name)
            quad_scaling_factors.append(scaling_factor)

        model_scaled.update()
        # Add scaling factors to scaled model
        model_scaled._quad_scaling_factors = quad_scaling_factors

    # Compute total time and capture final statistics
    total_time = time.time() - total_start_time
    final_stats = _capture_model_stats(model_scaled)

    # Store scaling time as model attribute
    model_scaled._scaling_time = total_time

    # Emit scaling footer with final stats
    logger.info(f"Scaling completed in {total_time:.2f} seconds")
    logger.info("Scaled Model Ranges:")
    logger.info(_extract_range_stats(final_stats))

    # Remove logging handlers added for this call
    for _h in _log_handlers:
        _h.close()
        logger.removeHandler(_h)

    # Restore original OutputFlag if we changed it
    if original_output_flag == 0 and needs_output:
        model.setParam("OutputFlag", 0)

    return model_scaled
