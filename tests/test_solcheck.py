import io
import pathlib
import contextlib
from unittest import TestCase

import gurobipy as gp
from gurobi_modelanalyzer import SolCheck

here = pathlib.Path(__file__).parent


class TestSolCheck(TestCase):
    def setUp(self):
        self.env = gp.Env()
        self.afiro_model = gp.read(str(here / "dataset" / "afiro.mps"), env=self.env)
        self.misc07_model = gp.read(str(here / "dataset" / "misc07.mps"), env=self.env)

    def tearDown(self):
        self.afiro_model.close()
        self.misc07_model.close()
        self.env.close()

    def test_suboptimal_solution(self):
        m = self.afiro_model
        sol = {m.getVarByName("X01"): 78, m.getVarByName("X22"): 495}
        sc = SolCheck(m)

        captured_output = io.StringIO()
        with contextlib.redirect_stdout(captured_output):
            sc.test_sol(sol)
        self.assertIn(
            "Solution is feasible for feasibility tolerance of 1e-06",
            captured_output.getvalue(),
        )

        sc.optimize()

        for v in sol.keys():
            self.assertLess(sol[v], v.X)
        self.assertGreater(abs(sc.obj_diff), 1e-6)

    # Test an infeasible solution
    def test_infeasible_solution(self):
        m = self.misc07_model
        sol = {m.getVarByName("COL260"): 2400.5}
        sc = SolCheck(m)

        captured_output = io.StringIO()
        with contextlib.redirect_stdout(captured_output):
            sc.test_sol(sol)

        self.assertIn(
            "Model is infeasible",
            captured_output.getvalue(),
        )
        self.assertEqual(sc.Status, gp.GRB.INFEASIBLE)

        # Fix violations

        captured_output = io.StringIO()
        with contextlib.redirect_stdout(captured_output):
            sc.inf_repair()

        self.assertIn(
            "Relaxing to find smallest violation from fixed solution",
            captured_output.getvalue(),
        )
        self.assertIn(
            "Fixed values are 1.5 from a feasible solution",
            captured_output.getvalue(),
        )

        abs_sum_violations = 0
        violation_rows = []
        violation_rows_names = []

        for c in m.getConstrs():
            if abs(c._Violation) > 1e-4:
                abs_sum_violations += abs(c._Violation)
                violation_rows.append(c._Violation)
                violation_rows_names.append(c.ConstrName)

        self.assertEqual(abs_sum_violations, 1.5)
        self.assertEqual(violation_rows, [-0.5, 1.0])
        self.assertEqual(violation_rows_names, ["ROW001", "ROW074"])
