import io
import os
import pathlib
import unittest
import warnings
from contextlib import redirect_stderr

import gurobipy as gp
from gurobipy import GRB
import numpy as np
import scipy.sparse

from gurobi_modelanalyzer import scale_model
from gurobi_modelanalyzer.scaling.methods import (
    ModelData,
    _threshold_small_coefficients,
    _compute_constraint_violations,
)
from gurobi_modelanalyzer.scaling.scaled_wrappers import (
    ScaledModel,
    ScaledVar,
    ScaledQConstr,
)

here = pathlib.Path(__file__).parent
cwd = pathlib.Path(os.getcwd())


# ── Shared helpers ────────────────────────────────────────────────────────


def _make_env():
    """Return a silent Gurobi environment."""
    env = gp.Env(params={"OutputFlag": 0})
    return env


def _tiny_lp(env):
    """
    A small LP whose coefficient matrix has very different magnitudes so
    that non-trivial scaling factors are produced.

      min    x + 1000 y
      s.t.   0.001 x + 1000 y >= 1
             x, y >= 0
    """
    m = gp.Model(env=env)
    x = m.addVar(lb=0, name="x")
    y = m.addVar(lb=0, name="y")
    m.addConstr(0.001 * x + 1000 * y >= 1, name="c1")
    m.setObjective(x + 1000 * y, GRB.MINIMIZE)
    m.update()
    return m


def _tiny_qp(env):
    """
    A small QP (quadratic objective, linear constraint).

      min   x^2 + y^2
      s.t.  x + y >= 1
            x, y >= 0
    """
    m = gp.Model(env=env)
    x = m.addVar(lb=0, name="x")
    y = m.addVar(lb=0, name="y")
    m.addConstr(x + y >= 1, name="c1")
    m.setObjective(x * x + y * y, GRB.MINIMIZE)
    m.update()
    return m


# ── scale_model API ───────────────────────────────────────────────────────


class TestScaleModelAPI(unittest.TestCase):
    """Tests for the scale_model() public API on a real LP."""

    def setUp(self):
        self.env = _make_env()
        self.model = gp.read(str(here / "dataset" / "afiro.mps"), env=self.env)
        self.model.setParam("OutputFlag", 0)

    def tearDown(self):
        self.model.close()
        self.env.close()

    # ── all three LP methods produce a ScaledModel ────────────────────────

    def test_equilibration_returns_scaled_model(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        self.assertIsInstance(ms, ScaledModel)
        ms.close()

    def test_geometric_mean_returns_scaled_model(self):
        ms = scale_model(
            self.model, "geometric_mean", env=self.env, scaling_log_to_console=0
        )
        self.assertIsInstance(ms, ScaledModel)
        ms.close()

    def test_arithmetic_mean_returns_scaled_model(self):
        ms = scale_model(
            self.model, "arithmetic_mean", env=self.env, scaling_log_to_console=0
        )
        self.assertIsInstance(ms, ScaledModel)
        ms.close()

    # ── structural checks ─────────────────────────────────────────────────

    def test_scaling_matrices_match_model_dimensions(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        n_vars = self.model.NumVars
        n_constrs = self.model.NumConstrs
        self.assertEqual(ms.ColScaling.shape, (n_vars, n_vars))
        self.assertEqual(ms.RowScaling.shape, (n_constrs, n_constrs))
        ms.close()

    def test_scaling_time_is_positive(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        self.assertIsNotNone(ms.ScalingTime)
        self.assertGreater(ms.ScalingTime, 0.0)
        ms.close()

    def test_variable_names_match_original(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        orig_names = [v.VarName for v in self.model.getVars()]
        scaled_names = [v.VarName for v in ms.getVars()]
        self.assertEqual(orig_names, scaled_names)
        ms.close()

    def test_constraint_names_match_original(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        orig_names = [c.ConstrName for c in self.model.getConstrs()]
        scaled_names = [c.ConstrName for c in ms.getConstrs()]
        self.assertEqual(orig_names, scaled_names)
        ms.close()

    def test_var_and_constr_counts_preserved(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        self.assertEqual(ms.NumVars, self.model.NumVars)
        self.assertEqual(ms.NumConstrs, self.model.NumConstrs)
        ms.close()

    # ── logging ───────────────────────────────────────────────────────────

    def test_silent_when_log_to_console_is_0(self):
        with redirect_stderr(io.StringIO()) as console:
            ms = scale_model(
                self.model, "equilibration", env=self.env, scaling_log_to_console=0
            )
            self.assertIsInstance(ms, ScaledModel)
            ms.close()
        output = console.getvalue()
        self.assertNotIn("Scaling Method:", output)

    def test_silent_when_log_to_console_is_1(self):
        with redirect_stderr(io.StringIO()) as console:
            ms = scale_model(
                self.model, "equilibration", env=self.env, scaling_log_to_console=1
            )
            self.assertIsInstance(ms, ScaledModel)
            ms.close()
        output = console.getvalue()
        self.assertIn("Scaling Method:", output)

    def test_log_written_to_file(self):
        log_path = cwd / "test_scaling_output.log"
        try:
            ms = scale_model(
                self.model,
                "equilibration",
                env=self.env,
                scaling_log=str(log_path),
                scaling_log_to_console=0,
            )
            self.assertTrue(log_path.exists())
            self.assertGreater(log_path.stat().st_size, 0)
            ms.close()
        finally:
            if log_path.exists():
                log_path.unlink()


# ── post-solve unscaling ──────────────────────────────────────────────────


class TestScalingFactorAttributes(unittest.TestCase):
    """scaling_factor is accessible on original var/constr objects and ScaledVar."""

    def setUp(self):
        self.env = _make_env()
        self.model = _tiny_lp(self.env)

    def tearDown(self):
        self.model.close()
        self.env.close()

    def _scale(self, **kwargs):
        return scale_model(
            self.model, "equilibration", scaling_log_to_console=0, **kwargs
        )

    def test_var_scaling_factor_matches_col_scaling_diagonal(self):
        ms = self._scale()
        ms.optimize()
        col_diag = ms.ColScaling.diagonal()
        for sv, expected in zip(ms.getVarsUnscaled(), col_diag):
            self.assertAlmostEqual(sv.scaling_factor, float(expected))
        ms.close()

    def test_constr_scaling_factor_matches_row_scaling_diagonal(self):
        ms = self._scale()
        row_diag = ms.RowScaling.diagonal()
        for sc, expected in zip(ms.getConstrsUnscaled(), row_diag):
            self.assertAlmostEqual(sc.scaling_factor, float(expected))
        ms.close()

    def test_all_var_scaling_factors_are_positive(self):
        ms = self._scale()
        ms.optimize()
        for sv in ms.getVarsUnscaled():
            self.assertGreater(sv.scaling_factor, 0.0)
        ms.close()

    def test_all_constr_scaling_factors_are_positive(self):
        ms = self._scale()
        for sc in ms.getConstrsUnscaled():
            self.assertGreater(sc.scaling_factor, 0.0)
        ms.close()

    def test_scaled_var_wrapper_scaling_factor_matches_original(self):
        ms = self._scale()
        ms.optimize()
        col_diag = ms.ColScaling.diagonal()
        for sv, expected in zip(ms.getVarsUnscaled(), col_diag):
            self.assertAlmostEqual(sv.scaling_factor, float(expected))
        ms.close()

    def test_qconstr_scaling_factor_is_set(self):
        """QP-constraint scaling factors are accessible via getQConstrsUnscaled()."""
        env = _make_env()
        try:
            m = _tiny_qcp(env)
            ms = scale_model(m, "equilibration", scaling_log_to_console=0)
            for sqc in ms.getQConstrsUnscaled():
                self.assertGreater(sqc.scaling_factor, 0.0)
            ms.close()
        finally:
            m.close()
            env.close()


class TestScaleModelSolve(unittest.TestCase):
    """Solve the scaled model and verify the unscaled solution is feasible."""

    def setUp(self):
        self.env = _make_env()
        self.model = gp.read(str(here / "dataset" / "afiro.mps"), env=self.env)
        self.model.setParam("OutputFlag", 0)

    def tearDown(self):
        self.model.close()
        self.env.close()

    def test_solve_yields_feasible_unscaled_solution(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()
        self.assertEqual(ms.Status, GRB.OPTIMAL)

        ms.computeUnscVio()
        self.assertIsNotNone(ms.MaxUnscVio)
        self.assertLess(ms.MaxUnscVio, 1e-4)
        ms.close()

    def test_getVarsUnscaled_returns_ScaledVar_objects(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()

        scaled_vars = ms.getVarsUnscaled()
        self.assertEqual(len(scaled_vars), self.model.NumVars)
        for sv in scaled_vars:
            self.assertIsInstance(sv, ScaledVar)
            # Xunsc must exist and be a finite number
            self.assertTrue(np.isfinite(sv.Xunsc))
        ms.close()


# ── _scale = 0 opt-out ────────────────────────────────────────────────────


class TestScaleOptOut(unittest.TestCase):
    """Verify that _scale=0 excludes variables/constraints from scaling."""

    def setUp(self):
        self.env = _make_env()
        self.model = _tiny_lp(self.env)

    def tearDown(self):
        self.model.close()
        self.env.close()

    def test_var_with_scale_0_not_column_scaled(self):
        """Column scaling factor for a var with _scale=0 must be exactly 1."""
        self.model.getVars()[0]._scale = 0
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        col_diag = ms.ColScaling.diagonal()
        self.assertAlmostEqual(col_diag[0], 1.0, places=10)
        ms.close()

    def test_constr_with_scale_0_not_row_scaled(self):
        """Row scaling factor for a constraint with _scale=0 must be exactly 1."""
        self.model.getConstrs()[0]._scale = 0
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        row_diag = ms.RowScaling.diagonal()
        self.assertAlmostEqual(row_diag[0], 1.0, places=10)
        ms.close()


# ── internal functions ────────────────────────────────────────────────────


class TestInternalFunctions(unittest.TestCase):
    """Unit tests for internal helpers in methods.py."""

    # _threshold_small_coefficients ----------------------------------------

    def test_threshold_zeros_small_values_in_sparse_matrix(self):
        data = scipy.sparse.csr_matrix(
            np.array([[1e-14, 1.0], [0.5, 1e-20]], dtype=float)
        )
        result = _threshold_small_coefficients(data)
        self.assertEqual(result[0, 0], 0.0)
        self.assertEqual(result[0, 1], 1.0)
        self.assertEqual(result[1, 0], 0.5)
        self.assertEqual(result[1, 1], 0.0)

    def test_threshold_zeros_small_values_in_dense_array(self):
        arr = np.array([1e-14, 0.5, 1e-20, 1.0])
        result = _threshold_small_coefficients(arr)
        self.assertEqual(result[0], 0.0)
        self.assertAlmostEqual(result[1], 0.5)
        self.assertEqual(result[2], 0.0)
        self.assertAlmostEqual(result[3], 1.0)

    def test_threshold_does_not_alter_large_values(self):
        arr = np.array([1.0, 100.0, 0.001])
        result = _threshold_small_coefficients(arr)
        np.testing.assert_array_almost_equal(result, arr)

    # _compute_constraint_violations ----------------------------------------

    def test_constraint_violation_for_infeasible_point(self):
        env = _make_env()
        m = gp.Model(env=env)
        x = m.addVar(lb=0, ub=10, name="x")
        m.addConstr(x >= 5, name="c_ge")  # violated when x=0
        m.addConstr(x <= 3, name="c_le")  # not violated when x=0
        m.setObjective(x, GRB.MINIMIZE)
        m.update()

        result = _compute_constraint_violations(m, [0.0])
        self.assertAlmostEqual(result["constraints"]["c_ge"], 5.0)
        self.assertAlmostEqual(result["constraints"]["c_le"], 0.0)
        m.close()
        env.close()

    def test_bound_violation_detected(self):
        env = _make_env()
        m = gp.Model(env=env)
        x = m.addVar(lb=2.0, ub=5.0, name="x")
        m.addConstr(x >= 0, name="c1")
        m.setObjective(x, GRB.MINIMIZE)
        m.update()

        # x=1 violates lb=2 → bound violation = 1
        result = _compute_constraint_violations(m, [1.0])
        self.assertAlmostEqual(result["bounds"]["x"], 1.0)
        m.close()
        env.close()

    # ModelData.from_gurobi_model -------------------------------------------

    def test_modeldata_captures_correct_dimensions(self):
        env = _make_env()
        m = gp.Model(env=env)
        x = m.addVar(lb=0, name="x")
        y = m.addVar(lb=0, name="y")
        m.addConstr(x + y >= 1, name="c1")
        m.setObjective(x + y, GRB.MINIMIZE)
        m.update()

        data = ModelData.from_gurobi_model(m)
        self.assertEqual(data.constr_matrix.shape, (1, 2))
        self.assertEqual(data.var_names, ["x", "y"])
        self.assertEqual(data.constr_names, ["c1"])
        np.testing.assert_array_equal(data.obj_vector, [1.0, 1.0])
        m.close()
        env.close()


def _tiny_qcp(env):
    """
    A small QCP: linear objective, one quadratic constraint.

      min   x + y
      s.t.  x^2 + y^2 <= 2
            x, y >= 0
    """
    m = gp.Model(env=env)
    x = m.addVar(lb=0, name="x")
    y = m.addVar(lb=0, name="y")
    m.addQConstr(x * x + y * y <= 2, name="qc1")
    m.setObjective(x + y, GRB.MINIMIZE)
    m.update()
    return m


# ── QP scaling ────────────────────────────────────────────────────────────


class TestQPScaling(unittest.TestCase):
    """Tests for models with a quadratic objective."""

    def setUp(self):
        self.env = _make_env()
        self.model = _tiny_qp(self.env)

    def tearDown(self):
        self.model.close()
        self.env.close()

    def test_non_equilibration_method_raises_warning(self):
        with self.assertWarns(UserWarning):
            ms = scale_model(
                self.model, "geometric_mean", env=self.env, scaling_log_to_console=0
            )
            ms.close()

    def test_equilibration_returns_scaled_model(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        self.assertIsInstance(ms, ScaledModel)
        ms.close()

    def test_qp_scaled_model_is_solvable(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()
        self.assertEqual(ms.Status, GRB.OPTIMAL)
        ms.close()


# ── QCP scaling ───────────────────────────────────────────────────────────


class TestQCPScaling(unittest.TestCase):
    """Tests for models with quadratic *constraints* (no quadratic objective)."""

    def setUp(self):
        self.env = _make_env()
        self.model = _tiny_qcp(self.env)

    def tearDown(self):
        self.model.close()
        self.env.close()

    # ── all three methods run on a QCP ────────────────────────────────────

    def test_equilibration_returns_scaled_model(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        self.assertIsInstance(ms, ScaledModel)
        ms.close()

    def test_geometric_mean_returns_scaled_model(self):
        ms = scale_model(
            self.model, "geometric_mean", env=self.env, scaling_log_to_console=0
        )
        self.assertIsInstance(ms, ScaledModel)
        ms.close()

    def test_arithmetic_mean_returns_scaled_model(self):
        ms = scale_model(
            self.model, "arithmetic_mean", env=self.env, scaling_log_to_console=0
        )
        self.assertIsInstance(ms, ScaledModel)
        ms.close()

    # ── qconstr count and suffix ──────────────────────────────────────────

    def test_qconstr_count_preserved(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        self.assertEqual(ms.NumQConstrs, self.model.NumQConstrs)
        ms.close()

    def test_quad_scaling_factors_stored(self):
        """_quad_scaling_factors must contain one entry per quadratic constraint."""
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        self.assertTrue(hasattr(ms, "_quad_scaling_factors"))
        self.assertEqual(len(ms._quad_scaling_factors), self.model.NumQConstrs)
        ms.close()

    def test_quad_scaling_factors_are_positive(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        for sf in ms._quad_scaling_factors:
            self.assertGreater(sf, 0.0)
        ms.close()

    # ── getQConstrsUnscaled ───────────────────────────────────────────────

    def test_getQConstrsUnscaled_returns_ScaledQConstr_objects(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()
        ms.computeUnscVio()

        qconstrs = ms.getQConstrsUnscaled()
        self.assertEqual(len(qconstrs), self.model.NumQConstrs)
        for qc in qconstrs:
            self.assertIsInstance(qc, ScaledQConstr)
        ms.close()

    # ── solve + unscaling ─────────────────────────────────────────────────

    def test_qcp_scaled_model_is_solvable(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()
        self.assertEqual(ms.Status, GRB.OPTIMAL)
        ms.close()

    def test_qcp_unscaled_solution_satisfies_quadratic_constraint(self):
        """
        After solving the scaled QCP, computeUnscVio should report zero
        violation for the quadratic constraint at the unscaled solution.
        """
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()
        self.assertEqual(ms.Status, GRB.OPTIMAL)

        ms.computeUnscVio()
        self.assertIsNotNone(ms.MaxUnscVio)
        self.assertLess(ms.MaxUnscVio, 1e-4)
        ms.close()

    def test_qcp_qconstr_unscaled_violation_is_populated(self):
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()
        ms.computeUnscVio()

        for qc in ms.getQConstrsUnscaled():
            self.assertIsNotNone(qc.UnscViolation)
            self.assertIsInstance(qc.UnscViolation, float)
        ms.close()

    # ── _scale = 0 opt-out for quadratic constraints ──────────────────────

    def test_qconstr_with_scale_0_has_scaling_factor_1(self):
        """A quadratic constraint with _scale=0 must get a scaling factor of 1.0."""
        self.model.getQConstrs()[0]._scale = 0
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        self.assertAlmostEqual(ms._quad_scaling_factors[0], 1.0, places=10)
        ms.close()


# ── init_scaling feature ──────────────────────────────────────────────────


class TestInitScaling(unittest.TestCase):
    """
    Tests for init_scaling=1 (user-provided only) and
    init_scaling=2 (warmstart).
    """

    def setUp(self):
        self.env = _make_env()
        # Use the ill-conditioned LP so non-trivial scaling factors emerge
        self.model = _tiny_lp(self.env)

    def tearDown(self):
        self.model.close()
        self.env.close()

    # ── invalid values ────────────────────────────────────────────────────

    def test_invalid_init_scaling_raises_value_error(self):
        for bad in (-1, 3, 99):
            with self.assertRaises(ValueError):
                scale_model(
                    self.model,
                    "equilibration",
                    init_scaling=bad,
                    scaling_log_to_console=0,
                )

    # ── mode 0 (default): _init_scaling is ignored ────────────────────────

    def test_mode0_ignores_init_scaling_attribute(self):
        """With init_scaling=0, _init_scaling attributes have no effect."""
        self.model.getVars()[0]._init_scaling = 42.0
        ms_plain = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        ms_attr = scale_model(
            self.model, "equilibration", init_scaling=0, scaling_log_to_console=0
        )
        col0_plain = ms_plain.ColScaling.diagonal()[0]
        col0_attr = ms_attr.ColScaling.diagonal()[0]
        self.assertAlmostEqual(col0_plain, col0_attr, places=10)
        ms_plain.close()
        ms_attr.close()

    # ── mode 1: use only _init_scaling, no algorithm ──────────────────────

    def test_mode1_uses_var_init_scaling_as_col_factor(self):
        """ColScaling diagonal must equal the _init_scaling values."""
        USER_FACTOR = 3.7
        for v in self.model.getVars():
            v._init_scaling = USER_FACTOR
        ms = scale_model(
            self.model, "equilibration", init_scaling=1, scaling_log_to_console=0
        )
        for factor in ms.ColScaling.diagonal():
            self.assertAlmostEqual(factor, USER_FACTOR, places=10)
        ms.close()

    def test_mode1_uses_constr_init_scaling_as_row_factor(self):
        """RowScaling diagonal must equal the _init_scaling values."""
        USER_FACTOR = 0.5
        for c in self.model.getConstrs():
            c._init_scaling = USER_FACTOR
        ms = scale_model(
            self.model, "equilibration", init_scaling=1, scaling_log_to_console=0
        )
        for factor in ms.RowScaling.diagonal():
            self.assertAlmostEqual(factor, USER_FACTOR, places=10)
        ms.close()

    def test_mode1_partial_init_scaling_defaults_to_1(self):
        """Variables without _init_scaling set should keep factor 1.0."""
        self.model.getVars()[0]._init_scaling = 5.0
        ms = scale_model(
            self.model, "equilibration", init_scaling=1, scaling_log_to_console=0
        )
        col_diag = ms.ColScaling.diagonal()
        self.assertAlmostEqual(col_diag[0], 5.0, places=10)
        self.assertAlmostEqual(col_diag[1], 1.0, places=10)
        ms.close()

    def test_mode1_init_scaling_overrides_scale_0(self):
        """_init_scaling takes priority over _scale=0."""
        USER_FACTOR = 2.0
        v = self.model.getVars()[0]
        v._scale = 0
        v._init_scaling = USER_FACTOR
        ms = scale_model(
            self.model, "equilibration", init_scaling=1, scaling_log_to_console=0
        )
        self.assertAlmostEqual(ms.ColScaling.diagonal()[0], USER_FACTOR, places=10)
        ms.close()

    def test_mode1_scaled_model_is_solvable(self):
        for v in self.model.getVars():
            v._init_scaling = 10.0
        ms = scale_model(
            self.model, "equilibration", init_scaling=1, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()
        self.assertEqual(ms.Status, GRB.OPTIMAL)
        ms.close()

    # ── mode 2: warmstart — algorithm runs on top of _init_scaling ────────

    def test_mode2_col_scaling_differs_from_mode0(self):
        """Warmstart with non-trivial init factors should change final scaling."""
        for v in self.model.getVars():
            v._init_scaling = 10.0
        ms0 = scale_model(
            self.model, "equilibration", init_scaling=0, scaling_log_to_console=0
        )
        ms2 = scale_model(
            self.model, "equilibration", init_scaling=2, scaling_log_to_console=0
        )
        # At least one column factor should differ
        differ = not np.allclose(ms0.ColScaling.diagonal(), ms2.ColScaling.diagonal())
        self.assertTrue(differ)
        ms0.close()
        ms2.close()

    def test_mode2_scale_0_with_init_locking_keeps_factor(self):
        """
        In mode 2, a variable with _scale=0 (excluded from the algorithm)
        AND _init_scaling=k should have final col factor exactly k.
        """
        USER_FACTOR = 4.0
        v = self.model.getVars()[0]
        v._scale = 0
        v._init_scaling = USER_FACTOR
        ms = scale_model(
            self.model, "equilibration", init_scaling=2, scaling_log_to_console=0
        )
        self.assertAlmostEqual(ms.ColScaling.diagonal()[0], USER_FACTOR, places=10)
        ms.close()

    def test_mode2_scaled_model_is_solvable(self):
        for v in self.model.getVars():
            v._init_scaling = 2.0
        ms = scale_model(
            self.model, "equilibration", init_scaling=2, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()
        self.assertEqual(ms.Status, GRB.OPTIMAL)
        ms.close()

    def test_mode2_unscaled_solution_is_feasible(self):
        for v in self.model.getVars():
            v._init_scaling = 2.0
        ms = scale_model(
            self.model, "equilibration", init_scaling=2, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()
        ms.computeUnscVio()
        self.assertLess(ms.MaxUnscVio, 1e-4)
        ms.close()


class TestInitScalingQCP(unittest.TestCase):
    """init_scaling tests specific to QCP quadratic constraints."""

    def setUp(self):
        self.env = _make_env()
        self.model = _tiny_qcp(self.env)

    def tearDown(self):
        self.model.close()
        self.env.close()

    def test_mode1_qconstr_init_scaling_used_as_row_factor(self):
        """In mode 1, qconstr._init_scaling must become the row scaling factor."""
        USER_FACTOR = 0.25
        self.model.getQConstrs()[0]._init_scaling = USER_FACTOR
        ms = scale_model(
            self.model, "equilibration", init_scaling=1, scaling_log_to_console=0
        )
        self.assertAlmostEqual(ms._quad_scaling_factors[0], USER_FACTOR, places=10)
        ms.close()

    def test_mode1_qconstr_init_scaling_overrides_scale_0(self):
        """_init_scaling overrides _scale=0 for QCP row factor in mode 1."""
        USER_FACTOR = 3.0
        qc = self.model.getQConstrs()[0]
        qc._scale = 0
        qc._init_scaling = USER_FACTOR
        ms = scale_model(
            self.model, "equilibration", init_scaling=1, scaling_log_to_console=0
        )
        self.assertAlmostEqual(ms._quad_scaling_factors[0], USER_FACTOR, places=10)
        ms.close()

    def test_mode2_qconstr_scale_0_with_init_locks_factor(self):
        """In mode 2, _scale=0 + _init_scaling locks the QCP row factor."""
        USER_FACTOR = 2.0
        qc = self.model.getQConstrs()[0]
        qc._scale = 0
        qc._init_scaling = USER_FACTOR
        ms = scale_model(
            self.model, "equilibration", init_scaling=2, scaling_log_to_console=0
        )
        self.assertAlmostEqual(ms._quad_scaling_factors[0], USER_FACTOR, places=10)
        ms.close()


# ── MIQCP: computeUnscVio and computeUnscObj do not raise TypeError ───────


def _tiny_miqcp(env):
    """
    A small MIQCP: integer variable, quadratic constraint, quadratic objective.

      min   x^2 + y^2 + z^2
      s.t.  x^2 + y^2 <= 4
            x + y + z >= 1
            x, y >= 0  (continuous)
            z in {0, 1} (binary)
    """
    m = gp.Model(env=env)
    x = m.addVar(lb=0, name="x")
    y = m.addVar(lb=0, name="y")
    z = m.addVar(vtype=GRB.BINARY, name="z")
    m.addQConstr(x * x + y * y <= 4, name="qc1")
    m.addConstr(x + y + z >= 1, name="c1")
    m.setObjective(x * x + y * y + z * z, GRB.MINIMIZE)
    m.update()
    return m


class TestMIQCPScaling(unittest.TestCase):
    """
    Regression tests for TypeError when float() is called on a numpy matrix
    returned by sparse @ dense quadratic expressions in computeUnscVio and
    computeUnscObj.
    """

    def setUp(self):
        self.env = _make_env()
        self.model = _tiny_miqcp(self.env)

    def tearDown(self):
        self.model.close()
        self.env.close()

    def test_compute_unsc_vio_does_not_raise(self):
        """computeUnscVio must not raise TypeError on a solved MIQCP."""
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()
        self.assertEqual(ms.Status, GRB.OPTIMAL)
        # Must not raise: "only 0-dimensional arrays can be converted to
        # Python scalars"
        ms.computeUnscVio()
        self.assertIsNotNone(ms.MaxUnscVio)
        self.assertIsInstance(ms.MaxUnscVio, float)
        ms.close()

    def test_compute_unsc_vio_value_is_near_zero(self):
        """Unscaled violations must be negligible at an optimal solution."""
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()
        self.assertEqual(ms.Status, GRB.OPTIMAL)
        ms.computeUnscVio()
        self.assertLess(ms.MaxUnscVio, 1e-4)
        ms.close()

    def test_compute_unsc_obj_does_not_raise(self):
        """computeUnscObj must not raise TypeError on a solved MIQCP."""
        ms = scale_model(
            self.model, "equilibration", env=self.env, scaling_log_to_console=0
        )
        ms.setParam("OutputFlag", 0)
        ms.optimize()
        self.assertEqual(ms.Status, GRB.OPTIMAL)
        ms.computeUnscObj()
        self.assertIsNotNone(ms.UnscObjVal)
        self.assertIsInstance(ms.UnscObjVal, float)
        ms.close()


if __name__ == "__main__":
    unittest.main()
