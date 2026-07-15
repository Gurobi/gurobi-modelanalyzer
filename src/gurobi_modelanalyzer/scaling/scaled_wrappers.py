import gurobipy as gp
from gurobipy import GRB
import scipy.sparse
import numpy as np
from typing import List, Union

from .methods import _compute_constraint_violations


class ScaledVar:
    """
    Wrapper around a Gurobi variable that provides access to unscaled
    values and bound violations.
    """

    def __init__(self, gurobi_var, col_scaling_factor):
        self._var = gurobi_var
        self._col_scaling_factor = col_scaling_factor
        self.UnscBoundViolation = None  # Will be set by computeUnscVio

    @property
    def X(self):
        """Scaled variable value"""
        return self._var.X

    @property
    def Xunsc(self):
        """Unscaled variable value: x = s * y"""
        return self._col_scaling_factor * self._var.X

    @property
    def scaling_factor(self):
        """Column scaling factor applied to this variable."""
        return self._col_scaling_factor

    def __getattr__(self, name):
        """Forward all other attributes to the underlying Gurobi variable"""
        return getattr(self._var, name)


class _ScaledConstraintBase:
    """
    Shared base for ScaledConstr and ScaledQConstr.
    Stores the wrapped Gurobi object, tracks unscaled violation, and
    forwards unknown attribute lookups to the wrapped object.
    """

    def __init__(self, wrapped):
        self._wrapped = wrapped
        self._unsc_violation = None

    @property
    def scaling_factor(self):
        """Row scaling factor applied to this constraint."""
        return self._wrapped._scaling_factor

    @property
    def UnscViolation(self):
        """Unscaled constraint violation"""
        return self._unsc_violation

    @UnscViolation.setter
    def UnscViolation(self, value):
        self._unsc_violation = value

    def __getattr__(self, name):
        """Forward all other attributes to the underlying Gurobi object"""
        return getattr(self._wrapped, name)


class ScaledConstr(_ScaledConstraintBase):
    """
    Wrapper around a Gurobi constraint that provides access to
    unscaled violations.
    """

    def __init__(self, gurobi_constr):
        super().__init__(gurobi_constr)


class ScaledQConstr(_ScaledConstraintBase):
    """
    Wrapper around a Gurobi quadratic constraint that provides
    access to unscaled violations.
    """

    def __init__(self, gurobi_qconstr):
        super().__init__(gurobi_qconstr)


class ScaledModel(gp.Model):
    """
    A Gurobi model with scaling information attached.
    Allows easy access to unscaled variable values after optimization.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._col_scaling = None
        self._row_scaling = None
        self._scaled_vars = None
        self._scaled_constrs = None
        self._scaled_qconstrs = None
        self._constraint_violations = None
        self._max_unsc_vio = None
        self._original_model = None

    def getVarsUnscaled(self):
        """
        Get list of ScaledVar objects that have both X and Xunsc attributes.

        Returns:
        --------
        List[ScaledVar]
            List of wrapped variables with unscaling capabilities
        """
        if self._scaled_vars is None:
            if self._col_scaling is None:
                raise ValueError("No scaling information available.")

            gurobi_vars = self.getVars()
            col_scaling_diag = self._col_scaling.diagonal()

            self._scaled_vars = [
                ScaledVar(var, col_scaling_diag[i]) for i, var in enumerate(gurobi_vars)
            ]

        return self._scaled_vars

    def _get_cached_wrappers(self, cache_attr, get_method, wrapper_cls):
        """Lazy-initialise and return a list of constraint wrapper objects."""
        if getattr(self, cache_attr) is None:
            if self._original_model is None:
                raise ValueError(
                    "Original model not stored. Cannot access constraint wrappers."
                )
            objects = get_method(self._original_model)
            setattr(self, cache_attr, [wrapper_cls(o) for o in objects])
        return getattr(self, cache_attr)

    def getConstrsUnscaled(self):
        """
        Get list of ScaledConstr objects with UnscViolation attributes.

        Returns:
        --------
        List[ScaledConstr]
            List of wrapped constraints with violation tracking
        """
        return self._get_cached_wrappers(
            "_scaled_constrs", lambda m: m.getConstrs(), ScaledConstr
        )

    def getQConstrsUnscaled(self):
        """
        Get list of ScaledQConstr objects with UnscViolation attributes.

        Returns:
        --------
        List[ScaledQConstr]
            List of wrapped quadratic constraints with violation tracking
        """
        return self._get_cached_wrappers(
            "_scaled_qconstrs", lambda m: m.getQConstrs(), ScaledQConstr
        )

    def computeUnscVio(self):
        """
        Compute unscaled constraint and bound violations using the
        unscaled variable values.
        Stores violations in constraint wrapper objects and variable
        wrapper objects,
        and tracks maximum violation.

        After calling this method, violations can be accessed via:
        - constraint.UnscViolation for each constraint
        - var.UnscBoundViolation for each variable's bound violation
        - model_scaled.MaxUnscVio for the maximum violation across
          all constraints and bounds
        - model_scaled.MaxUnscConstrVio for max constraint violation only
        - model_scaled.MaxUnscBoundVio for max bound violation only
        """
        if self._original_model is None:
            raise ValueError(
                "Original model not set. ScaledModel must be created via scale_model()."
            )

        # Get unscaled solution values
        unscaled_vars = self.getVarsUnscaled()
        unscaled_values = [var.Xunsc for var in unscaled_vars]

        # Compute constraint and bound violations using the original model
        violations = _compute_constraint_violations(
            self._original_model, unscaled_values
        )

        # Store violations in dictionaries
        self._constraint_violations = violations["constraints"]
        self._bound_violations = violations["bounds"]

        # Store violations in constraint wrappers
        scaled_constrs = self.getConstrsUnscaled()
        for scaled_constr in scaled_constrs:
            constr_name = scaled_constr.ConstrName
            scaled_constr.UnscViolation = self._constraint_violations.get(
                constr_name, 0.0
            )

        # Store violations in quadratic constraint wrappers
        if self._original_model.NumQConstrs > 0:
            scaled_qconstrs = self.getQConstrsUnscaled()
            for scaled_qconstr in scaled_qconstrs:
                qconstr_name = scaled_qconstr.QCName
                scaled_qconstr.UnscViolation = self._constraint_violations.get(
                    qconstr_name, 0.0
                )

        # Store bound violations in variable wrappers
        for i, var in enumerate(unscaled_vars):
            var_name = var.VarName
            var.UnscBoundViolation = self._bound_violations.get(var_name, 0.0)

        # Compute and store maximum violations
        all_constraint_vios = list(self._constraint_violations.values())
        all_bound_vios = list(self._bound_violations.values())

        self._max_unsc_constr_vio = (
            max(all_constraint_vios) if all_constraint_vios else 0.0
        )
        self._max_unsc_bound_vio = max(all_bound_vios) if all_bound_vios else 0.0
        self._max_unsc_vio = max(self._max_unsc_constr_vio, self._max_unsc_bound_vio)

    @property
    def MaxUnscVio(self):
        """
        Get the maximum unscaled violation across all constraints and bounds.

        Returns:
        --------
        float
            Maximum violation, or None if not computed
        """
        return getattr(self, "_max_unsc_vio", None)

    @property
    def MaxUnscConstrVio(self):
        """
        Get the maximum unscaled constraint violation (linear and
        quadratic constraints only).

        Returns:
        --------
        float
            Maximum constraint violation, or None if not computed
        """
        return getattr(self, "_max_unsc_constr_vio", None)

    @property
    def MaxUnscBoundVio(self):
        """
        Get the maximum unscaled bound violation.

        Returns:
        --------
        float
            Maximum bound violation, or None if not computed
        """
        return getattr(self, "_max_unsc_bound_vio", None)

    @property
    def ScalingTime(self):
        """
        Get the time taken to scale the model.

        Returns:
        --------
        float
            Scaling time in seconds, or None if not available
        """
        return getattr(self, "_scaling_time", None)

    @property
    def ColScaling(self):
        """
        Get the column scaling matrix.

        Returns:
        --------
        scipy.sparse.csr_matrix
            Diagonal matrix with column scaling factors, or None if
            not available
        """
        return getattr(self, "_col_scaling", None)

    @property
    def RowScaling(self):
        """
        Get the row scaling matrix.

        Returns:
        --------
        scipy.sparse.csr_matrix
            Diagonal matrix with row scaling factors, or None if not available
        """
        return getattr(self, "_row_scaling", None)

    def computeUnscObj(self):
        """
        Compute the unscaled objective value using original model coefficients
        and unscaled variable values.

        This method should be called after optimization to get the objective
        value in the original (unscaled) space. The result is stored in
        :py:attr:`UnscObjVal`.
        """
        if self._original_model is None:
            raise ValueError(
                "Original model not set. ScaledModel must be created via scale_model()."
            )

        # Get unscaled solution values
        unscaled_vars = self.getVarsUnscaled()
        unscaled_values = np.array([var.Xunsc for var in unscaled_vars])

        # Get original objective coefficients
        orig_obj = np.array(self._original_model.getAttr("Obj"))

        # Compute linear objective contribution
        linear_obj = np.dot(orig_obj, unscaled_values)

        # Check for quadratic objective
        q_matrix = self._original_model.getQ()
        if q_matrix.nnz > 0:
            # Quadratic contribution: x^T q x
            # q is upper triangular, need full symmetric form
            q_full = q_matrix + q_matrix.T - scipy.sparse.diags(q_matrix.diagonal())
            quad_obj = float(
                np.asarray(unscaled_values @ q_full @ unscaled_values).flat[0]
            )
        else:
            quad_obj = 0.0

        self._unsc_obj_val = linear_obj + quad_obj

    @property
    def UnscObjVal(self):
        """
        Get the unscaled objective value.

        This is the objective value computed using original model coefficients
        and unscaled variable values.

        Returns:
        --------
        float
            Unscaled objective value, or None if not computed.
            Call :py:meth:`computeUnscObj` first.
        """
        return getattr(self, "_unsc_obj_val", None)

    @property
    def OriginalModel(self):
        """
        The original (unscaled) Gurobi model that was passed to
        :func:`scale_model`.

        Returns:
        --------
        gp.Model or None
        """
        return self._original_model

    def write_scaling(self, path: str, lock_factors: bool = True) -> None:
        """
        Export the scaling factors to a ``.scl`` file.

        The file can be passed back to ``gurobi_cls`` via ``--scaling-file``
        or to :func:`scale_model` after parsing with
        :func:`~gurobi_modelanalyzer.scaling.cli._parse_scaling_file`, to
        reproduce or continue from the same scaling.

        All variables and constraints are written, including those with a
        factor of 1.0. When ``lock_factors=True``, a factor of 1.0 with
        ``lock_flag=0`` explicitly locks that object at the identity scaling
        and prevents the algorithm from modifying it on re-import.

        Parameters
        ----------
        path : str
            Output file path (conventionally with ``.scl`` extension).
        lock_factors : bool, optional
            If ``True`` (default), all entries are written with
            ``lock_flag = 0``, meaning the factors are kept fixed when the
            file is re-imported and no further refinement is performed.
            If ``False``, ``lock_flag = 1`` is written, so the factors act
            as a warmstart and the scaling algorithm may adjust them further.
        """
        if self._col_scaling is None or self._row_scaling is None:
            raise ValueError(
                "No scaling information available. "
                "ScaledModel must be created via scale_model()."
            )
        if self._original_model is None:
            raise ValueError(
                "Original model reference not available. "
                "ScaledModel must be created via scale_model()."
            )

        lock_flag = 0 if lock_factors else 1
        col_diag = self._col_scaling.diagonal()
        row_diag = self._row_scaling.diagonal()

        orig_vars = self._original_model.getVars()
        orig_constrs = self._original_model.getConstrs()

        quad_factors = getattr(self, "_quad_scaling_factors", None)
        has_qconstrs = bool(quad_factors)

        with open(path, "w") as fh:
            fh.write("# Gurobi scaling output file\n")
            fh.write("GRB_SCL_FILE_VERSION 1\n")
            fh.write("\n")

            fh.write("SECTION VARS\n")
            for var, factor in zip(orig_vars, col_diag):
                fh.write(f"{var.VarName}  {factor:.15g}  {lock_flag}\n")

            fh.write("\n")
            fh.write("SECTION CONSTRS\n")
            for constr, factor in zip(orig_constrs, row_diag):
                fh.write(f"{constr.ConstrName}  {factor:.15g}  {lock_flag}\n")

            if has_qconstrs:
                fh.write("\n")
                fh.write("SECTION QCONSTRS\n")
                orig_qconstrs = self._original_model.getQConstrs()
                for qconstr, factor in zip(orig_qconstrs, quad_factors):
                    fh.write(
                        f"{qconstr.QCName}  {factor:.15g}  {lock_flag}\n"
                    )
