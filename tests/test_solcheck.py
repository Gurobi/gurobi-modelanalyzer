import io
import os
import pathlib
import hashlib
import subprocess
import contextlib
from unittest import TestCase

import gurobipy as gp
from gurobi_modelanalyzer import SolCheck

here = pathlib.Path(__file__).parent
cwd = pathlib.Path(os.getcwd())


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

        for c in m.getConstrs():
            if abs(c._Violation) > 1e-4:
                abs_sum_violations += abs(c._Violation)

        self.assertEqual(abs_sum_violations, 1.5)


class TestSolCheckCLI(TestCase):
    def setUp(self):
        self.afirofix_md5 = "f265f0c9554f2b480092757a83e47207"

    def tearDown(self):
        output_files = [
            cwd.joinpath("afirofix.sol"),
            cwd.joinpath("misc07fix.vio"),
            cwd.joinpath("misc07fix.sol"),
            cwd.joinpath("misc07fix.ilp"),
        ]
        for file in output_files:
            if file.exists():
                os.remove(file)

    def test_suboptimal_solution(self):
        afirofix_file = cwd.joinpath("afirofix.sol")
        self.assertFalse(afirofix_file.exists())
        cmd = (
            f"gurobi_solcheck --model {str(here / 'dataset' / 'afiro.mps')} "
            + f"--sol {str(here / 'dataset' / 'afiro.sol')} --result afirofix"
        )
        result = subprocess.getoutput(cmd)
        self.assertIn(
            "Solution is feasible for feasibility tolerance of 1e-06",
            result,
        )
        self.assertIn("Difference: -5.0613", result)
        self.assertTrue(afirofix_file.exists())

        md5_ret = None
        with open(afirofix_file, "rb") as f:
            d = f.read()
            md5_ret = hashlib.md5(d).hexdigest()

        # broken on Windows, disabled for now
        # self.assertEqual(self.afirofix_md5, md5_ret)

    def test_suboptimal_json_solution(self):
        afirofix_file = cwd.joinpath("afirofix.sol")
        self.assertFalse(afirofix_file.exists())
        cmd = (
            f"gurobi_solcheck --model {str(here / 'dataset' / 'afiro.mps')} "
            + f"--sol {str(here / 'dataset' / 'afiro.json')} --result afirofix"
        )
        result = subprocess.getoutput(cmd)
        self.assertIn(
            "Solution is feasible for feasibility tolerance of 1e-06",
            result,
        )
        self.assertIn("Difference: -5.0613", result)
        self.assertTrue(afirofix_file.exists())

        md5_ret = None
        with open(afirofix_file, "rb") as f:
            d = f.read()
            md5_ret = hashlib.md5(d).hexdigest()

        # broken on Windows, disabled for now
        # self.assertEqual(self.afirofix_md5, md5_ret)

    def test_infeasible_solution(self):
        misc07fix_vio_file = cwd.joinpath("misc07fix.vio")
        self.assertFalse(misc07fix_vio_file.exists())
        cmd = (
            f"gurobi_solcheck --model {str(here / 'dataset' / 'misc07.mps')} "
            + f"--sol {str(here / 'dataset' / 'misc07.sol')} --result misc07fix"
        )
        result = subprocess.getoutput(cmd)
        self.assertTrue(misc07fix_vio_file.exists())

        # Check log for expected output
        self.assertIn(
            "Model is infeasible",
            result,
        )
        self.assertIn(
            "Solution is infeasible for feasibility tolerance of 1e-06", result
        )
        self.assertIn(
            "Relaxing to find smallest violation from fixed solution",
            result,
        )
        self.assertIn("Fixed values are 1.5 from a feasible solution", result)

    def test_infeasible_solution_infmethodV(self):
        misc07fix_sol_file = cwd.joinpath("misc07fix.sol")
        self.assertFalse(misc07fix_sol_file.exists())
        cmd = (
            f"gurobi_solcheck --model {str(here / 'dataset' / 'misc07.mps')} "
            + f"--sol {str(here / 'dataset' / 'misc07.sol')} --result misc07fix "
            + "--infmethod V "
        )
        result = subprocess.getoutput(cmd)
        self.assertTrue(misc07fix_sol_file.exists())

        # Check log for expected output
        self.assertIn("Model is infeasible", result)
        self.assertIn(
            "Relaxing to find smallest violation from fixed solution",
            result,
        )
        self.assertIn("Fixed values are 409.5 from a feasible solution", result)

    def test_infeasible_solution_infmethodI(self):
        misc07fix_ilp_file = cwd.joinpath("misc07fix.ilp")
        self.assertFalse(misc07fix_ilp_file.exists())
        cmd = (
            f"gurobi_solcheck --model {str(here / 'dataset' / 'misc07.mps')} "
            + f"--sol {str(here / 'dataset' / 'misc07.sol')} --result misc07fix "
            + "--infmethod I "
        )
        result = subprocess.getoutput(cmd)
        self.assertTrue(misc07fix_ilp_file.exists())

        # Check log for expected output
        self.assertIn("Model is infeasible", result)
        self.assertIn(
            "Solution is infeasible for feasibility tolerance of 1e-06", result
        )
        self.assertIn(
            "Computing Irreducible Inconsistent Subsystem (IIS)",
            result,
        )
        self.assertIn("IIS computed: 1 constraints, 2 bounds", result)
