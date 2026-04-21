import gurobipy as gp
from gurobipy import GRB
import scipy.sparse
import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Union
import time
import logging

logger = logging.getLogger("gurobi_modelanalyzer.scaling")


def _extract_range_stats(stats: str) -> str:
    """
    Extract only range-related lines from model statistics.

    Parameters:
    -----------
    stats : str
        Full model statistics string from printStats()

    Returns:
    --------
    str
        Only the lines containing range information
    """
    range_lines = [line for line in stats.split("\n") if "range" in line.lower()]
    return "\n".join(range_lines)


def _print_scaling_log(
    method: str,
    original_stats: str,
    scale_passes: int,
    iteration_logs: List[Dict],
    total_time: float,
    final_stats: str,
    scaling_time_limit: float = float("inf"),
    scale_conv_tol: float = 1e-4,
    mode: str = "final",
):
    """
    Emit formatted scaling log entries via the module logger.

    Attach handlers to ``logging.getLogger('gurobi_modelanalyzer.scaling')``
    to control where output is written.  ``scale_model`` manages handler
    setup and teardown automatically based on its ``scaling_log`` and
    ``scaling_log_to_console`` parameters.

    Parameters:
    -----------
    method : str
        Scaling method name
    original_stats : str
        Statistics output from original model's printStats()
    scale_passes : int
        Number of scale passes requested
    iteration_logs : List[Dict]
        List of dictionaries containing iteration information
    total_time : float
        Total elapsed time for scaling
    final_stats : str
        Statistics output from scaled model's printStats()
    scaling_time_limit : float, optional
        Time limit for scaling (default: inf)
    mode : str, optional
        'header'    = emit header and iteration-table header
        'iteration' = emit a single iteration row
        'final'     = emit complete log as one block (default)
    """
    if mode == "header":
        logger.info(f"Scaling Method: {method}")
        logger.info(f"Scale Passes:   {scale_passes}")
        logger.info(f"Conv. Tol.:     {scale_conv_tol:.6e}")
        if scaling_time_limit != float("inf"):
            logger.info(f"Time Limit:     {scaling_time_limit:.2f} seconds")
        logger.info("Original Model Statistics:")
        logger.info(original_stats.rstrip())
        logger.info(f"{'Scale Pass':<12} {'Max Factor Dev.':<17} {'Time (s)':<15}")
        return

    if mode == "iteration":
        if len(iteration_logs) > 0:
            log = iteration_logs[-1]
            pass_num = log.get("pass", "-")
            rel_change = log.get("rel_change", 0.0)
            iter_time = log.get("time", 0.0)
            if isinstance(rel_change, float):
                logger.info(f"{pass_num:<12} {rel_change:<17.6e} {iter_time:<15.2f}")
            else:
                logger.info(f"{pass_num:<12} {rel_change:<17} {iter_time:<15.2f}")
        return

    # mode == 'final': emit complete log as a single block
    log_lines = []
    log_lines.append(f"Scaling Method: {method}")
    log_lines.append(f"Scale Passes:   {scale_passes}")
    log_lines.append(f"Conv. Tol.:     {scale_conv_tol:.6e}")
    if scaling_time_limit != float("inf"):
        log_lines.append(f"Time Limit:     {scaling_time_limit:.2f} seconds")
    log_lines.append("Original Model Statistics:")
    log_lines.append(original_stats.rstrip())
    log_lines.append(f"{'Scale Pass':<12} {'Max Factor Dev.':<17} {'Time (s)':<15}")
    for log in iteration_logs:
        pass_num = log.get("pass", "-")
        rel_change = log.get("rel_change", 0.0)
        iter_time = log.get("time", 0.0)
        if isinstance(rel_change, float):
            log_lines.append(f"{pass_num:<12} {rel_change:<17.6e} {iter_time:<15.2f}")
        else:
            log_lines.append(f"{pass_num:<12} {rel_change:<17} {iter_time:<15.2f}")
    log_lines.append(f"Scaling completed in {total_time:.2f} seconds")
    log_lines.append("Scaled Model Ranges:")
    log_lines.append(_extract_range_stats(final_stats))
    logger.info("\n".join(log_lines))


@dataclass
class ModelData:
    """
    Encapsulates Gurobi model data extracted prior to scaling.
    """

    constr_matrix: scipy.sparse.csr_matrix = field(default=None)
    rhs_vector: np.ndarray = field(default=None)
    constr_sense: List[str] = field(default=None)
    obj_vector: np.ndarray = field(default=None)
    ub_vector: np.ndarray = field(default=None)
    lb_vector: np.ndarray = field(default=None)
    var_types: List[str] = field(default=None)
    var_names: List[str] = field(default=None)
    constr_names: List[str] = field(default=None)

    @classmethod
    def from_gurobi_model(cls, model):
        """
        Create ModelData from a Gurobi model.

        Parameters:
        -----------
        model : gp.Model
            Gurobi model to extract data from

        Returns:
        --------
        ModelData
            ModelData object containing all model information
        """

        return cls(
            constr_matrix=model.getA(),
            rhs_vector=np.array(model.getAttr("RHS")),
            constr_sense=model.getAttr("Sense"),
            obj_vector=np.array(model.getAttr("Obj")),
            ub_vector=np.array(model.getAttr("UB")),
            lb_vector=np.array(model.getAttr("LB")),
            var_types=model.getAttr("VType"),
            var_names=[var.VarName for var in model.getVars()],
            constr_names=[constr.ConstrName for constr in model.getConstrs()],
        )


def _row_scale_factor(row_data: np.ndarray, method: str) -> float:
    """Compute a row scaling factor for the given method."""
    if method == "equilibration":
        return 1.0 / np.max(row_data)
    elif method == "arithmetic_mean":
        return 1.0 / np.mean(row_data)
    else:  # geometric_mean
        return 1.0 / np.sqrt(np.min(row_data) * np.max(row_data))


def _col_scale_factor(col_data: np.ndarray, method: str) -> float:
    """Compute a column scaling factor for the given method."""
    if method == "equilibration":
        return 1.0 / np.max(col_data)
    elif method == "arithmetic_mean":
        return 1.0 / np.mean(col_data)
    else:  # geometric_mean
        return 1.0 / np.sqrt(np.min(col_data) * np.max(col_data))


def _iterative_scaling(
    constr_matrix: scipy.sparse.csr_matrix,
    cols_to_scale: List[int],
    rows_to_scale: List[int],
    scale_passes: int,
    scale_conv_tol: float,
    method: str,
    scaling_time_limit: float = float("inf"),
) -> Tuple[
    scipy.sparse.csr_matrix,
    scipy.sparse.csr_matrix,
    scipy.sparse.csr_matrix,
    List[Dict],
]:
    """
    Shared iterative row/column scaling loop used by equilibration,
    geometric_mean, and arithmetic_mean.

    Parameters:
    -----------
    constr_matrix : scipy.sparse.csr_matrix
        The constraint matrix to scale
    cols_to_scale : List[int]
        Column indices to scale
    rows_to_scale : List[int]
        Row indices to scale
    scale_passes : int
        Maximum number of scaling iterations
    scale_conv_tol : float
        Relative tolerance for convergence check
    method : str
        Scaling method ('equilibration', 'geometric_mean', or
        'arithmetic_mean')
    scaling_time_limit : float, optional
        Time limit in seconds (default: inf - no limit)

    Returns:
    --------
    Tuple[
        scipy.sparse.csr_matrix, scipy.sparse.csr_matrix,
        scipy.sparse.csr_matrix, List[Dict]]
        - scaled_constr_matrix: The scaled constraint matrix
        - row_scaling_total: Cumulative row scaling matrix (diagonal)
        - col_scaling_total: Cumulative column scaling matrix (diagonal)
        - iteration_logs: List of iteration information dictionaries
    """
    scaled_constr_matrix = constr_matrix
    num_rows, num_cols = constr_matrix.shape
    row_scaling_total = scipy.sparse.eye(num_rows, format="csr")
    col_scaling_total = scipy.sparse.eye(num_cols, format="csr")
    iteration_logs = []
    total_elapsed_time = 0.0

    rows_to_scale_set = set(rows_to_scale)
    cols_to_scale_set = set(cols_to_scale)

    scaled_constr_matrix_csc = scaled_constr_matrix.tocsc()

    for completed_scale_passes in range(scale_passes):
        iter_start_time = time.time()
        scaled_constr_matrix = scaled_constr_matrix_csc.tocsr()

        # Compute row scaling factors (skip excluded rows)
        row_factors = np.ones(num_rows)
        for i in rows_to_scale_set:
            row_data = np.abs(scaled_constr_matrix.getrow(i).data)
            if len(row_data) > 0:
                row_factors[i] = _row_scale_factor(row_data, method)

        row_scaling_iter = scipy.sparse.diags(row_factors)
        scaled_constr_matrix = row_scaling_iter @ scaled_constr_matrix
        row_scaling_total = row_scaling_iter @ row_scaling_total

        scaled_constr_matrix_csc = scaled_constr_matrix.tocsc()

        # Compute column scaling factors via direct CSC array access
        col_factors_full = np.ones(num_cols)
        csc_data = np.abs(scaled_constr_matrix_csc.data)
        csc_indptr = scaled_constr_matrix_csc.indptr
        for j in range(num_cols):
            if j in cols_to_scale_set:
                start_idx = csc_indptr[j]
                end_idx = csc_indptr[j + 1]
                if end_idx > start_idx:
                    col_factors_full[j] = _col_scale_factor(
                        csc_data[start_idx:end_idx], method
                    )

        col_scaling_iter = scipy.sparse.diags(col_factors_full)
        scaled_constr_matrix_csc = scaled_constr_matrix_csc @ col_scaling_iter
        col_scaling_total = col_scaling_total @ col_scaling_iter

        rel_change = max(
            np.max(np.abs(row_factors - 1.0)) if row_factors.size > 0 else 0.0,
            np.max(np.abs(col_factors_full - 1.0))
            if col_factors_full.size > 0
            else 0.0,
        )

        iter_time = time.time() - iter_start_time
        total_elapsed_time += iter_time
        iteration_logs.append(
            {
                "pass": completed_scale_passes + 1,
                "rel_change": rel_change,
                "time": total_elapsed_time,
            }
        )
        _print_scaling_log("", "", 0, iteration_logs, 0.0, "", mode="iteration")

        if rel_change < scale_conv_tol:
            break
        if total_elapsed_time >= scaling_time_limit:
            break

    scaled_constr_matrix = scaled_constr_matrix_csc.tocsr()
    return (scaled_constr_matrix, row_scaling_total, col_scaling_total, iteration_logs)


def _scale_single_qconstr(
    qconstr: gp.QConstr,
    model: gp.Model,
    col_scaling: scipy.sparse.csr_matrix,
    skip_row_scale: bool = False,
    row_scale_override: float = None,
    scaling_lb: float = 1e-8,
) -> Tuple[scipy.sparse.csr_matrix, np.ndarray, str, float, float, str]:
    """
    Compute scaling for a single quadratic constraint.

    Parameters:
    -----------
    qconstr : gp.QConstr
        Quadratic constraint to scale
    model : gp.Model
        Original model containing the constraint
    col_scaling : scipy.sparse.csr_matrix
        Column scaling matrix
    skip_row_scale : bool, optional
        If True, skip the row-level normalisation (scaling_factor=1).
        Column scaling is still applied. Use when constr._scale=0
        (default: False)
    row_scale_override : float, optional
        When set, use this value directly as the row scaling factor,
        overriding both skip_row_scale and the norm-based computation.
        Used to honour constr._init_scaling (default: None)
    scaling_lb : float, optional
        Lower bound for the row scaling factor to avoid extreme values
        (default: 1e-8)

    Returns:
    --------
    Tuple containing:
        - qc_scaled: Scaled quadratic matrix
        - q_scaled: Scaled linear vector (as flattened numpy array)
        - sense: Constraint sense
        - rhs_scaled: Scaled RHS
        - scaling_factor: Computed scaling factor
        - name: Constraint name
    """
    # Extract data from constraint
    qc, q = model.getQCMatrices(qconstr)
    qc = col_scaling @ qc @ col_scaling
    q = col_scaling @ q
    # Convert q to dense array and flatten (getQCMatrices returns sparse
    # column matrix)
    if scipy.sparse.issparse(q):
        q = np.asarray(q.todense()).flatten()
    else:
        q = np.asarray(q).flatten()
    rhs = qconstr.QCRHS
    sense = qconstr.QCSense
    name = qconstr.QCName + "_scaled"

    # Compute Frobenius norm efficiently for upper triangular matrix
    # For symmetric matrix: ||A||_F^2 = sum(A_ij^2) for all i,j
    # For upper triangular q: ||q + q^T - diag(q)||_F^2 =
    # 2*sum(Q_ij^2) - sum(Q_ii^2)
    # Faster: sqrt(2 * ||q||_F^2 - ||diag(q)||_2^2)
    if qc.nnz > 0:
        qc_norm_sq = np.sum(qc.data**2)  # ||q||_F^2
        qc_diag_norm_sq = np.sum(qc.diagonal() ** 2)  # ||diag(q)||_2^2
        qc_full_norm = np.sqrt(2.0 * qc_norm_sq - qc_diag_norm_sq)
    else:
        qc_full_norm = 0.0

    # Compute scaling factor for constraint
    q_norm = np.linalg.norm(q) if q.size > 0 else 0.0
    if row_scale_override is not None:
        scaling_factor = row_scale_override
    elif skip_row_scale:
        scaling_factor = 1.0
    else:
        scaling_factor = max(1.0 / max(qc_full_norm, q_norm, abs(rhs), 1.0), scaling_lb)

    qc_scaled = scaling_factor * qc
    q_scaled = scaling_factor * q
    rhs_scaled = scaling_factor * rhs

    return qc_scaled, q_scaled, sense, rhs_scaled, scaling_factor, name


def _threshold_small_coefficients(
    data: Union[scipy.sparse.spmatrix, np.ndarray],
    value_threshold: float = 1e-13,
) -> Union[scipy.sparse.csr_matrix, np.ndarray]:
    """
    Set coefficients below threshold to zero.
    Works for both sparse matrices and numpy arrays.

    Parameters:
    -----------
    data : scipy.sparse matrix or np.ndarray
        The data to threshold
    value_threshold : float, optional
        Absolute value threshold below which coefficients are set to
        zero (default: 1e-13)

    Returns:
    --------
    scipy.sparse.csr_matrix or np.ndarray
        Thresholded data in the same format as input
    """
    if scipy.sparse.issparse(data):
        A = data.tocsr().copy()
        mask = np.abs(A.data) >= value_threshold
        A.data = A.data * mask
        A.eliminate_zeros()
        return A
    else:
        # For numpy arrays
        result = data.copy()
        result[np.abs(result) < value_threshold] = 0.0
        return result


def equilibration(
    constr_matrix: scipy.sparse.csr_matrix,
    cols_to_scale: List[int],
    rows_to_scale: List[int],
    scale_passes: int,
    scale_conv_tol: float,
    scaling_time_limit: float = float("inf"),
) -> Tuple[
    scipy.sparse.csr_matrix,
    scipy.sparse.csr_matrix,
    scipy.sparse.csr_matrix,
    List[Dict],
]:
    """
    Scale constraint matrix using equilibration method.

    Applies iterative row and column scaling using mean-based equilibration.
    Row scaling uses the mean of absolute row values, column scaling uses
    the maximum of absolute column values.

    Parameters:
    -----------
    constr_matrix : scipy.sparse.csr_matrix
        The constraint matrix to scale
    cols_to_scale : List[int]
        List of column indices to scale (typically continuous variables)
    rows_to_scale : List[int]
        List of row indices to scale (constraints not in this list
        keep a row scaling factor of 1)
    scale_passes : int
        Maximum number of scaling iterations
    scale_conv_tol : float
        Relative tolerance for convergence check
    scaling_time_limit : float, optional
        Time limit in seconds for scaling iterations (default: inf - no limit)

    Returns:
    --------
    Tuple[
        scipy.sparse.csr_matrix, scipy.sparse.csr_matrix,
        scipy.sparse.csr_matrix, List[Dict]]
        - scaled_constr_matrix: The scaled constraint matrix
        - row_scaling_total: Cumulative row scaling matrix (diagonal)
        - col_scaling_total: Cumulative column scaling matrix (diagonal)
        - iteration_logs: List of iteration information dictionaries
    """
    return _iterative_scaling(
        constr_matrix,
        cols_to_scale,
        rows_to_scale,
        scale_passes,
        scale_conv_tol,
        "equilibration",
        scaling_time_limit,
    )


def geometric_mean(
    constr_matrix: scipy.sparse.csr_matrix,
    cols_to_scale: List[int],
    rows_to_scale: List[int],
    scale_passes: int,
    scale_conv_tol: float,
    scaling_time_limit: float = float("inf"),
) -> Tuple[
    scipy.sparse.csr_matrix,
    scipy.sparse.csr_matrix,
    scipy.sparse.csr_matrix,
    List[Dict],
]:
    """
    Scale constraint matrix using geometric mean method.

    Applies iterative row and column scaling using geometric mean of
    minimum and maximum absolute values in each row/column.

    Parameters:
    -----------
    constr_matrix : scipy.sparse.csr_matrix
        The constraint matrix to scale
    cols_to_scale : List[int]
        List of column indices to scale (typically continuous variables)
    rows_to_scale : List[int]
        List of row indices to scale (constraints not in this list
        keep a row scaling factor of 1)
    scale_passes : int
        Maximum number of scaling iterations
    scale_conv_tol : float
        Relative tolerance for convergence check
    scaling_time_limit : float, optional
        Time limit in seconds for scaling iterations (default: inf - no limit)

    Returns:
    --------
    Tuple[
        scipy.sparse.csr_matrix, scipy.sparse.csr_matrix,
        scipy.sparse.csr_matrix, List[Dict]]
        - scaled_constr_matrix: The scaled constraint matrix
        - row_scaling_total: Cumulative row scaling matrix (diagonal)
        - col_scaling_total: Cumulative column scaling matrix (diagonal)
        - iteration_logs: List of iteration information dictionaries
    """
    return _iterative_scaling(
        constr_matrix,
        cols_to_scale,
        rows_to_scale,
        scale_passes,
        scale_conv_tol,
        "geometric_mean",
        scaling_time_limit,
    )


def arithmetic_mean(
    constr_matrix: scipy.sparse.csr_matrix,
    cols_to_scale: List[int],
    rows_to_scale: List[int],
    scale_passes: int,
    scale_conv_tol: float,
    scaling_time_limit: float = float("inf"),
) -> Tuple[
    scipy.sparse.csr_matrix,
    scipy.sparse.csr_matrix,
    scipy.sparse.csr_matrix,
    List[Dict],
]:
    """
    Scale constraint matrix using arithmetic mean method.

    Applies iterative row and column scaling using arithmetic mean of
    absolute values in each row/column.

    Parameters:
    -----------
    constr_matrix : scipy.sparse.csr_matrix
        The constraint matrix to scale
    cols_to_scale : List[int]
        List of column indices to scale (typically continuous variables)
    rows_to_scale : List[int]
        List of row indices to scale (constraints not in this list
        keep a row scaling factor of 1)
    scale_passes : int
        Maximum number of scaling iterations
    scale_conv_tol : float
        Relative tolerance for convergence check
    scaling_time_limit : float, optional
        Time limit in seconds for scaling iterations (default: inf - no limit)

    Returns:
    --------
    Tuple[
        scipy.sparse.csr_matrix, scipy.sparse.csr_matrix,
        scipy.sparse.csr_matrix, List[Dict]]
        - scaled_constr_matrix: The scaled constraint matrix
        - row_scaling_total: Cumulative row scaling matrix (diagonal)
        - col_scaling_total: Cumulative column scaling matrix (diagonal)
        - iteration_logs: List of iteration information dictionaries
    """
    return _iterative_scaling(
        constr_matrix,
        cols_to_scale,
        rows_to_scale,
        scale_passes,
        scale_conv_tol,
        "arithmetic_mean",
        scaling_time_limit,
    )


def quad_equilibration(
    constr_matrix: scipy.sparse.csr_matrix,
    obj_vector: np.ndarray,
    q_matrix: scipy.sparse.coo_matrix,
    cols_to_scale: List[int],
    rows_to_scale: List[int],
    scale_passes: int,
    scale_conv_tol: float,
    scaling_lb: float = 1e-8,
    scaling_ub: float = 1e8,
    scaling_time_limit: float = float("inf"),
) -> Tuple[
    scipy.sparse.csr_matrix,
    scipy.sparse.csr_matrix,
    np.ndarray,
    scipy.sparse.csr_matrix,
    scipy.sparse.csr_matrix,
    List[Dict],
]:
    """
    Scale quadratic program using KKT-based equilibration method.

    Scales the constraint matrix, quadratic objective matrix (q), and linear
    objective vector jointly by building a KKT matrix and
    applying equilibration
    to both the variable/constraint scaling and the objective scaling.

    Parameters:
    -----------
    constr_matrix : scipy.sparse.csr_matrix
        The constraint matrix to scale
    obj_vector : np.ndarray
        Linear objective coefficient vector
    q_matrix : scipy.sparse.coo_matrix
        Quadratic objective matrix (Hessian)
    cols_to_scale : List[int]
        List of column indices to scale (typically continuous variables)
    rows_to_scale : List[int]
        List of row indices to scale (constraints not in this list
        keep a row scaling factor of 1)
    scale_passes : int
        Maximum number of scaling iterations
    scale_conv_tol : float
        Relative tolerance for convergence check
    scaling_lb : float, optional
        Lower bound for scaling factors (default: 1e-8)
    scaling_ub : float, optional
        Upper bound for scaling factors (default: 1e8)
    scaling_time_limit : float, optional
        Time limit in seconds for scaling iterations (default: inf - no limit)

    Returns:
    --------
    Tuple[
        scipy.sparse.csr_matrix, scipy.sparse.csr_matrix,
        np.ndarray, scipy.sparse.csr_matrix,
        scipy.sparse.csr_matrix, List[Dict]]
        - scaled_constr_matrix: The scaled constraint matrix
        - scaled_q_matrix: The scaled quadratic objective matrix
        - scaled_obj_vector: The scaled linear objective vector
        - row_scaling_total: Cumulative row scaling matrix (diagonal)
        - col_scaling_total: Cumulative column scaling matrix (diagonal)
        - iteration_logs: List of iteration information dictionaries
    """
    scaled_constr_matrix = constr_matrix.copy()
    scaled_q_matrix = q_matrix.copy()
    scaled_obj_vector = obj_vector.copy()

    # Initialize scaling matrices
    num_rows, num_cols = constr_matrix.shape
    diagonal_scaling_total = scipy.sparse.eye(num_rows + num_cols, format="csr")
    obj_scaling_factor_total = 1.0
    zero_block = scipy.sparse.csr_matrix((num_rows, num_rows))
    iteration_logs = []
    total_elapsed_time = 0.0

    for completed_scale_passes in range(scale_passes):
        iter_start_time = time.time()

        # Build KKT matrix from CURRENT scaled matrices and convert to CSC for
        # column operations
        kkt_matrix = scipy.sparse.bmat(
            [
                [scaled_q_matrix, scaled_constr_matrix.T],
                [scaled_constr_matrix, zero_block],
            ]
        ).tocsc()  # CSC format for efficient column access

        # Compute diagonal scaling factors using CSC direct access (much faster
        # than getcol)
        diagonal_factors = np.ones(num_rows + num_cols)
        cols_to_scale_set = set(cols_to_scale)
        rows_to_scale_set = set(rows_to_scale)

        # Access CSC internal arrays directly for speed
        kkt_data = np.abs(kkt_matrix.data)
        kkt_indptr = kkt_matrix.indptr

        for i in range(num_rows + num_cols):
            # Skip integer/binary variable columns
            if i < num_cols and i not in cols_to_scale_set:
                continue
            # Skip excluded constraint rows
            if i >= num_cols and (i - num_cols) not in rows_to_scale_set:
                continue

            # Get column data using direct array access
            start_idx = kkt_indptr[i]
            end_idx = kkt_indptr[i + 1]

            if end_idx > start_idx:  # Column has data
                col_data = kkt_data[start_idx:end_idx]
                max_val = np.max(col_data)
                scaling_factor = 1.0 / np.sqrt(max_val)
                diagonal_factors[i] = np.clip(scaling_factor, scaling_lb, scaling_ub)
        diagonal_scaling_iter = scipy.sparse.diags(diagonal_factors)

        # Extract column and row scaling from THIS iteration
        col_scaling_iter = scipy.sparse.diags(
            diagonal_scaling_iter.diagonal()[:num_cols]
        )
        row_scaling_iter = scipy.sparse.diags(
            diagonal_scaling_iter.diagonal()[num_cols:]
        )

        # Apply M equilibration scaling
        scaled_constr_matrix = (
            row_scaling_iter @ scaled_constr_matrix @ col_scaling_iter
        )
        scaled_q_matrix = col_scaling_iter @ scaled_q_matrix @ col_scaling_iter
        scaled_obj_vector = col_scaling_iter @ scaled_obj_vector

        # Compute cost scaling factor γ
        q_col_norms = []
        for j in range(num_cols):
            col_data = np.abs(scaled_q_matrix.getcol(j).data)
            if len(col_data) > 0:
                q_col_norms.append(np.max(col_data))

        denominator = max(
            np.mean(q_col_norms) if q_col_norms else 1.0,
            np.max(np.abs(scaled_obj_vector)) if scaled_obj_vector.size > 0 else 1.0,
            1.0,
        )
        obj_scaling_factor = 1.0 / denominator
        obj_scaling_factor = np.clip(obj_scaling_factor, scaling_lb, scaling_ub)

        # Apply cost scaling
        scaled_q_matrix = obj_scaling_factor * scaled_q_matrix
        scaled_obj_vector = obj_scaling_factor * scaled_obj_vector

        # Accumulate total scaling
        diagonal_scaling_total = diagonal_scaling_iter @ diagonal_scaling_total
        obj_scaling_factor_total *= obj_scaling_factor

        # Check convergence: max deviation of diagonal factors from 1
        # (same criterion as the LP case, matching Algorithm 2 in the paper)
        rel_change = np.max(np.abs(diagonal_factors - 1.0))

        iter_time = time.time() - iter_start_time
        total_elapsed_time += iter_time
        iteration_logs.append(
            {
                "pass": completed_scale_passes + 1,
                "rel_change": rel_change,
                "time": total_elapsed_time,
            }
        )

        # Emit iteration progress
        _print_scaling_log("", "", 0, iteration_logs, 0.0, "", mode="iteration")

        if rel_change < scale_conv_tol:
            break

        # Check time limit
        if total_elapsed_time >= scaling_time_limit:
            break

    # Extract final column and row scaling
    col_scaling_total = scipy.sparse.diags(diagonal_scaling_total.diagonal()[:num_cols])
    row_scaling_total = scipy.sparse.diags(diagonal_scaling_total.diagonal()[num_cols:])

    return (
        scaled_constr_matrix,
        scaled_q_matrix,
        scaled_obj_vector,
        row_scaling_total,
        col_scaling_total,
        iteration_logs,
    )


def _compute_violation(sense: str, lhs_value: float, rhs: float) -> float:
    """Compute constraint violation given sense, LHS value, and RHS."""
    if sense == GRB.LESS_EQUAL:
        return max(0.0, lhs_value - rhs)
    elif sense == GRB.GREATER_EQUAL:
        return max(0.0, rhs - lhs_value)
    elif sense == GRB.EQUAL:
        return abs(lhs_value - rhs)
    return 0.0


def _compute_constraint_violations(
    model: gp.Model,
    unscaled_variables: Union[List[float], np.ndarray],
) -> Dict[str, Dict[str, float]]:
    """
    Compute constraint and bound violations for a given candidate solution.

    Evaluates how much a candidate solution violates the constraints
    and variable
    bounds of the model. This is useful for checking solution quality or
    assessing infeasibility.

    Parameters:
    -----------
    model : gp.Model
        The Gurobi model containing constraints
    unscaled_variables : list or np.ndarray
        Solution values in the same order as model.getVars()

    Returns:
    --------
    Dict[str, Dict[str, float]]
        Dictionary with two keys:
        - 'constraints': dict mapping constraint names to violation values
        - 'bounds': dict mapping variable names to bound violation values

    Notes:
    ------
    Violation definitions:
        - For <= constraints: max(0, LHS - RHS)
        - For >= constraints: max(0, RHS - LHS)
        - For == constraints: |LHS - RHS|
        - For bounds: max(0, x - UB) + max(0, LB - x)
    """
    constraint_violations = {}
    bound_violations = {}
    vars_list = model.getVars()
    solution_dict = {var: val for var, val in zip(vars_list, unscaled_variables)}

    # Process bound violations
    for var, val in zip(vars_list, unscaled_variables):
        lb_vio = max(0.0, var.LB - val) if var.LB > -GRB.INFINITY else 0.0
        ub_vio = max(0.0, val - var.UB) if var.UB < GRB.INFINITY else 0.0
        bound_violations[var.VarName] = lb_vio + ub_vio  # One of the two will be zero

    # Process linear constraints
    for constr in model.getConstrs():
        # Get constraint properties
        sense = constr.Sense
        rhs = constr.RHS
        name = constr.ConstrName

        # Compute LHS value
        row = model.getRow(constr)
        lhs_value = 0.0
        for i in range(row.size()):
            var = row.getVar(i)
            coeff = row.getCoeff(i)
            lhs_value += coeff * solution_dict.get(var, 0.0)

        constraint_violations[name] = _compute_violation(sense, lhs_value, rhs)

    # Process quadratic constraints
    try:
        for qconstr in model.getQConstrs():
            # Get constraint properties
            sense = qconstr.QCSense
            rhs = qconstr.QCRHS
            name = qconstr.QCName

            # Get quadratic and linear parts
            q, c = model.getQCMatrices(qconstr)

            # Build solution vector in model variable order
            x = np.array([solution_dict.get(v, 0.0) for v in vars_list])

            # Compute LHS: x^T q x + c^T x
            # q is upper triangular, so we need to account for symmetry
            q_full = q + q.T - np.diag(q.diagonal())  # Make q symmetric
            lhs_value = float(np.asarray(x.T @ q_full @ x + c.T @ x).flat[0])

            constraint_violations[name] = _compute_violation(sense, lhs_value, rhs)

    except ImportError:
        # numpy not available, skip quadratic constraints
        num_qconstrs = model.NumQConstrs
        if num_qconstrs > 0:
            print(
                f"Warning: {num_qconstrs} quadratic constraint(s) "
                f"found but numpy is not available."
            )
            print("         Quadratic constraint violations will not be computed.")

    return {"constraints": constraint_violations, "bounds": bound_violations}
