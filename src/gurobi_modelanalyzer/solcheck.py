"""Use Gurobi to diagnose a given solution for feasibility and optimality"""

from gurobipy import GRB
import re
import gzip
import sys
import json
import itertools


def openc(fn, mode="r"):
    """Open a file with optional gzip compression"""
    return gzip.open(fn, mode + "t") if ext_check(fn, ("gz",)) else open(fn, mode)


def add_dots(fns):
    """Prepend dot to each filename in a collection (list)"""
    return map(lambda x: "." + x, fns)


def ext_check(fn, exts, compress=[]):
    """Check that a filename ends with a given list of extensions, with optional
    extensions for compression"""
    compress = [""] + list(add_dots(compress))
    return fn.lower().endswith(
        tuple(map("".join, itertools.product(add_dots(exts), compress)))
    )


class SolCheck:
    def __init__(self, model):
        self.model = model
        self.inf_method = None
        self.Status = GRB.LOADED
        self.obj_diff = None
        self.__solpat = re.compile("(\S+)\s+([+\-]?\d+\.?\d*([eE][+\-]?\d+)?)")

    def message(self, lines, pad=True):
        """Write lines using model logging"""
        if not isinstance(lines, list):
            # Embed string in a list
            lines = [lines]
        if pad:
            # Add blank lines before and after message
            lines.insert(0, "")
            lines.append("")
        for line in lines:
            self.model.message(line)

    def sol_fix(self, sol):
        """Fix variable bounds to solution values, and save original bounds"""
        self.sol = sol
        for var, val in self.sol.items():
            var._LB, var._UB = var.LB, var.UB
            var.LB, var.UB = val, val

    def sol_unfix(self):
        """Restore original variable bounds"""
        for var in self.model.getVars():
            try:
                var.LB, var.UB = var._LB, var._UB
            except AttributeError:
                pass

    def read_sol(self, fn, useZeros=False):
        """Read solution file"""
        if ext_check(fn, ("sol",), ("gz",)):
            # Handle SOL format
            sol = {}
            with openc(fn) as fp:
                for line in fp:
                    mat = self.__solpat.match(line)
                    if mat:
                        vname = mat.group(1)
                        val = float(mat.group(2))
                        sol[vname] = val
        elif ext_check(fn, ("json",), ("gz",)):
            # Handle JSON format
            with openc(fn) as fp:
                data = json.load(fp)
            sol = {v["VarName"]: v["X"] for v in data["Vars"]}
        else:
            raise RuntimeError(
                f"{fn} solution file must be either .sol(.gz) or .json(.gz)"
            )

        outsol = {}
        for vname, val in sol.items():
            if val != 0.0 or useZeros:
                var = self.model.getVarByName(vname)
                if var:
                    outsol[var] = val
                else:
                    self.message(f"Warning: {vname} not found in model", pad=False)

        return outsol

    def make_zero_sol(self):
        """Make zero solution"""
        return {v: 0.0 for v in self.model.getVars()}

    def write_vio_file(self, fn):
        """Write solution violation file"""
        with openc(fn, "w") as fp:
            fp.write(f"# Total constraint violation: {self.model.ObjVal}\n")
            for c in self.model.getConstrs():
                fp.write(f"{c.ConstrName} {c._Violation}\n")

    def write_result(self, fn):
        """Write output file"""
        if self.Status == GRB.LOADED:
            raise RuntimeError("Must test solution first")

        if self.Status == GRB.INFEASIBLE and self.inf_method == "I":
            if not ext_check(fn, ("ilp",), ("gz", "bz2", "7z", "zip")):
                fn += ".ilp"
        elif self.Status == GRB.INFEASIBLE and self.inf_method == "C":
            if not ext_check(fn, ("vio",), ("gz",)):
                fn += ".vio"
        else:
            if not ext_check(fn, ("sol",), ("gz", "bz2", "7z", "zip")):
                fn += ".sol"
        self.message([f"Writing result file {fn}"])
        if self.Status == GRB.INFEASIBLE and self.inf_method == "C":
            self.write_vio_file(fn)
        else:
            self.model.write(fn)

    def test_sol(self, sol):
        """Test a solution"""
        self.sol_fix(sol)
        self.model.optimize()
        self.Status = self.model.Status
        if self.Status == GRB.OPTIMAL:
            self.Status = GRB.SUBOPTIMAL
            feasMsg = "feasible"
        else:
            feasMsg = "infeasible"

        self.message(
            f"Solution is {feasMsg} for feasibility tolerance of {self.model.params.FeasibilityTol}"
        )

    def inf_explain(self):
        """Explain infeasible solution"""
        if self.Status == GRB.LOADED:
            raise RuntimeError("Must test solution first")

        # Force the solution into the IIS
        for var in self.sol.keys():
            var.IISLBForce = 1
            var.IISUBForce = 1

        self.model.computeIIS()
        self.inf_method = "I"

    def inf_repair(self, repairMethod="C", makeCopy=False):
        """Repair a solution or model"""
        if self.Status == GRB.LOADED:
            raise RuntimeError("Must test solution first")

        self.message("Relaxing to find smallest violation from fixed solution")
        origNumVars = self.model.NumVars

        if makeCopy:
            relax = self.model.copy()
            self.model._relax = relax
        else:
            relax = self.model

        if repairMethod == "C":  # relax constraints
            relax.feasRelax(
                0, False, None, None, None, relax.getConstrs(), relax.NumConstrs * [1]
            )
        elif repairMethod == "V":  # relax variables
            relax.feasRelax(
                0,
                False,
                list(self.sol.keys()),
                len(self.sol) * [1],
                len(self.sol) * [1],
                None,
                None,
            )
        else:
            raise RuntimeError(f"repairMethod {repairMethod} not recognized")

        relax.optimize()
        self.inf_method = repairMethod

        if repairMethod == "C":  # add attribute for constraint violations
            rval = iter(relax.getAttr("X", relax.getVars()[origNumVars:]))
            for c in self.model.getConstrs():
                c._Violation = 0.0
                if c.Sense in ("=", ">"):  # ArtP_
                    c._Violation += next(rval)
                if c.Sense in ("=", "<"):  # ArtN_
                    c._Violation -= next(rval)

        self.message(f"Fixed values are {relax.ObjVal} from a feasible solution")

    def optimize(self):
        """Optimize from test solution"""
        if self.Status == GRB.LOADED:
            raise RuntimeError("Must test solution first")

        self.message("Comparing quality with original solution")
        self.model._FixObjVal = self.model.ObjVal

        self.sol_unfix()
        self.model.optimize()

        self.obj_diff = self.model.ObjVal - self.model._FixObjVal

        # We use MIP gap even if the model is not a MIP to measure the sameness
        # of the objective
        if abs(self.obj_diff) < self.model.Params.MIPGap:
            self.Status = GRB.OPTIMAL
        else:
            self.Status = GRB.SUBOPTIMAL

        self.message(
            [
                "Objectives:",
                f"Fixed:      {self.model._FixObjVal:.4f}",
                f"Optimal:    {self.model.ObjVal:.4f}",
                f"Difference: {self.obj_diff:.4f}",
            ]
        )


# Command line code
def main_cli():
    import argparse
    import gurobipy as gp
    import questionary as qy

    parser = argparse.ArgumentParser(
        prog="gurobi_solcheck", description="Check if a solution is feasible"
    )

    parser.add_argument("--model", "-m", help="Model file")
    parser.add_argument("--sol", "-s", help="Solution file to test")
    parser.add_argument(
        "--result", "-r", help="Filename where the result should be written"
    )
    parser.add_argument(
        "--testonly", action="store_true", help="Test without solution diagnosis"
    )
    parser.add_argument(
        "--infmethod",
        choices=("C", "V", "I"),
        type=str.upper,
        default="C",
        help="Method to diagnose infeasible solution (C: Repair by adjusting RHS of constraints, V: Repair by adjusting variable values, I: Explain via IIS)",
    )
    parser.add_argument(
        "--usezeros",
        "-z",
        action="store_true",
        help="Use zero values in the solution file",
    )

    args = parser.parse_args()

    if args.model:
        mfn = args.model
    else:
        print("Entering interactive mode")
        mfn = qy.path(
            "Select a model file",
            validate=lambda fn: ext_check(
                fn, ("mps", "lp", "rew", "rlp"), ("gz", "bz2", "7z", "zip")
            ),
        ).ask()

    # Read model
    sc = SolCheck(gp.read(mfn))

    if args.model:
        sfn = args.sol
        useZeros = args.usezeros
    else:
        if qy.confirm("Use a solution file?").ask():
            sfn = qy.path(
                "Select a solution file",
                validate=lambda fn: ext_check(fn, ("sol", "json"), ("gz",)),
            ).ask()
            useZeros = not qy.confirm("Ignore zeros in solution file?").ask()
        else:
            sfn = ""

    if sfn:
        sol = sc.read_sol(sfn, useZeros)
    else:
        if args.model:
            sc.message("No solution file specified, using zero solution")
        sol = sc.make_zero_sol()

    # Check solution
    sc.test_sol(sol)

    if args.model:
        diagnose = not args.testonly
    else:
        diagnose = qy.confirm("Diagnose solution?").ask()

    if diagnose:
        if sc.Status == GRB.INFEASIBLE:
            if args.model:
                inf_method = args.infmethod
            else:
                inf_method = qy.select(
                    "Select method to diagnose infeasible solution",
                    choices=[
                        qy.Choice("Repair by adjusting RHS of constraints", "C"),
                        qy.Choice("Repair by adjusting variable values", "V"),
                        qy.Choice("Explain via IIS", "I"),
                    ],
                ).ask()

            if inf_method == "I":
                sc.inf_explain()
            else:
                sc.inf_repair(inf_method)
        else:
            sc.optimize()

        if args.model:
            rfn = args.result
        else:
            if qy.confirm("Write result file?").ask():
                rfn = qy.path("Select a result file").ask()
            else:
                rfn = ""

        if rfn:
            sc.write_result(rfn)

    sys.exit(sc.Status)
