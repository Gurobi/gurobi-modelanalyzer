"""Command-line interface for scaling and solving a Gurobi model."""

import argparse
import inspect
import io
import os
import sys
import tempfile


def _coerce_param_value(val: str):
    """Coerce a CLI string value to int or float if possible, else keep str."""
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        pass
    return val


def _capture_print_quality(model):
    """Return model.printQuality() output as a string."""
    old_stdout = sys.stdout
    sys.stdout = buf = io.StringIO()
    try:
        model.printQuality()
        return buf.getvalue()
    finally:
        sys.stdout = old_stdout


def _format_quality_stats(scaled_model):
    """Return a formatted quality-stats block covering scaled and unscaled."""
    lines = []

    # --- Scaled model ---
    lines.append("Solution quality - scaled model")
    lines.append(_capture_print_quality(scaled_model).rstrip())
    lines.append(f"  Objective value         : {scaled_model.ObjVal:.10g}")

    # --- Unscaled model ---
    scaled_model.computeUnscObj()
    scaled_model.computeUnscVio()
    lines.append("")
    lines.append("Solution quality - unscaled model (original variable space)")
    lines.append(f"  Objective value         : {scaled_model.UnscObjVal:.10g}")
    lines.append(
        f"  Max constraint violation: {scaled_model.MaxUnscConstrVio:.6e}"
    )
    lines.append(
        f"  Max bound violation     : {scaled_model.MaxUnscBoundVio:.6e}"
    )
    lines.append(
        f"  Max total violation     : {scaled_model.MaxUnscVio:.6e}"
    )
    lines.append("")

    return "\n".join(lines) + "\n"


def main_cli():
    import gurobipy as gp
    from gurobipy import GRB
    from gurobi_modelanalyzer.scaling import scale_model

    # Pull defaults directly from scale_model() so CLI always stays in sync
    _sig = inspect.signature(scale_model)
    _defaults = {
        k: v.default
        for k, v in _sig.parameters.items()
        if v.default is not inspect.Parameter.empty
    }

    parser = argparse.ArgumentParser(
        prog="gurobi_cls",
        description=(
            "Scale a Gurobi model and solve it.\n\n"
            "Reads the model file, applies scaling, solves the scaled model,\n"
            "and writes a unified log containing the scaling log, the Gurobi\n"
            "solve log, and solution quality statistics.\n\n"
            "Gurobi solver parameters can be passed as positional 'Param=Value'\n"
            "arguments before the model file, exactly like gurobi_cl:\n\n"
            "  gurobi_cls --method equilibration TimeLimit=60 Presolve=2 model.mps\n\n"
            "The log file defaults to 'gurobi.log'. Override with LogFile=path."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ---- Scaling-specific options ----
    parser.add_argument(
        "--method",
        "-m",
        choices=("equilibration", "geometric_mean", "arithmetic_mean"),
        default="equilibration",
        help=(
            "Scaling method: 'equilibration' (default), 'geometric_mean', or "
            "'arithmetic_mean'. Only 'equilibration' supports QP models."
        ),
    )
    parser.add_argument(
        "--scale-passes",
        type=int,
        default=_defaults["scale_passes"],
        metavar="N",
        help=f"Maximum number of scaling iterations (default: {_defaults['scale_passes']}).",
    )
    parser.add_argument(
        "--scale-conv-tol",
        type=float,
        default=_defaults["scale_conv_tol"],
        metavar="TOL",
        help=(
            f"Convergence tolerance: stop early when the max deviation of "
            f"scaling factors from 1 falls below this value "
            f"(default: {_defaults['scale_conv_tol']})."
        ),
    )
    parser.add_argument(
        "--scaling-lb",
        type=float,
        default=_defaults["scaling_lb"],
        metavar="LB",
        help=f"Lower bound for scaling factors (default: {_defaults['scaling_lb']}).",
    )
    parser.add_argument(
        "--scaling-ub",
        type=float,
        default=_defaults["scaling_ub"],
        metavar="UB",
        help=f"Upper bound for scaling factors (default: {_defaults['scaling_ub']}).",
    )
    parser.add_argument(
        "--value-threshold",
        type=float,
        default=_defaults["value_threshold"],
        metavar="THRESH",
        help=f"Coefficients below this threshold are set to zero (default: {_defaults['value_threshold']}).",
    )
    parser.add_argument(
        "--scaling-time-limit",
        type=float,
        default=_defaults["scaling_time_limit"],
        metavar="SEC",
        help="Time limit in seconds for scaling iterations (default: no limit).",
    )
    parser.add_argument(
        "--no-console-log",
        action="store_true",
        help="Suppress scaling log output to the console.",
    )
    parser.add_argument(
        "--scaling-file",
        metavar="PATH",
        default=None,
        help=(
            "Path to a .scl scaling input file specifying per-variable and "
            "per-constraint initial scaling factors. When provided, the "
            "scaling algorithm is run as a warmstart on top of the supplied "
            "factors (init_scaling mode 2). File format: section-based text "
            "with SECTION VARS / SECTION CONSTRS / SECTION QCONSTRS headers; "
            "each data line is 'name  factor  lock_flag' where lock_flag is "
            "0 (keep factor fixed) or 1 (use as warmstart, allow refinement)."
        ),
    )

    # ---- Positional: [Param=Value ...] model ----
    parser.add_argument(
        "args",
        nargs="+",
        metavar="[Param=Value ...] model",
        help=(
            "Optional Gurobi parameters in 'Param=Value' format "
            "(e.g. TimeLimit=60 Presolve=2), followed by the model file path."
        ),
    )

    parsed = parser.parse_args()

    # Split positionals: last = model file, rest = Gurobi Param=Value pairs
    positional = parsed.args
    model_path = positional[-1]
    gurobi_param_strs = positional[:-1]

    # Parse Param=Value strings
    gurobi_params = {}
    for p in gurobi_param_strs:
        if "=" not in p:
            parser.error(
                f"Gurobi parameters must be in 'Param=Value' format, got: {p!r}"
            )
        key, _, val = p.partition("=")
        gurobi_params[key] = val

    # Determine log file: respect LogFile= Gurobi param, default gurobi.log
    logfile_key = next((k for k in gurobi_params if k.lower() == "logfile"), None)
    logfile = gurobi_params[logfile_key] if logfile_key else "gurobi.log"

    # Read model
    try:
        model = gp.read(model_path)
    except Exception as exc:
        print(f"Error reading model file '{model_path}': {exc}", file=sys.stderr)
        sys.exit(1)

    # Apply scaling input file, if given
    from gurobi_modelanalyzer.scaling import read_scaling_file

    init_scaling_mode = 0
    if parsed.scaling_file is not None:
        try:
            read_scaling_file(parsed.scaling_file, model)
        except OSError as exc:
            print(
                f"Error reading scaling file '{parsed.scaling_file}': {exc}",
                file=sys.stderr,
            )
            sys.exit(1)
        init_scaling_mode = 2

    scaling_log_to_console = 0 if parsed.no_console_log else 1

    # Write the scaling log to a temp file so we can prepend it to the
    # final logfile after Gurobi's optimize() creates/overwrites it.
    tmp_fd, temp_scaling_log = tempfile.mkstemp(suffix=".log", prefix="gurobi_cls_")
    os.close(tmp_fd)

    opt_status = None
    quality_text = ""

    try:
        # Scale the model
        scaled_model = scale_model(
            model=model,
            method=parsed.method,
            scale_passes=parsed.scale_passes,
            scale_conv_tol=parsed.scale_conv_tol,
            scaling_lb=parsed.scaling_lb,
            scaling_ub=parsed.scaling_ub,
            value_threshold=parsed.value_threshold,
            scaling_time_limit=parsed.scaling_time_limit,
            scaling_log=temp_scaling_log,
            scaling_log_to_console=scaling_log_to_console,
            init_scaling=init_scaling_mode,
        )

        # Write scaling output file (locked factors, so the file can be
        # used to reproduce the same scaling exactly on re-import).
        model_stem = os.path.splitext(os.path.basename(model_path))[0]
        scl_output_path = model_stem + ".scl"
        scaled_model.write_scaling(scl_output_path, lock_factors=True)
        print(f"Scaling factors written to: {scl_output_path}")

        # Apply Gurobi solver parameters to the scaled model.
        # Skip LogFile here; we set it explicitly below so Gurobi writes
        # its solve log to our target logfile during optimize().
        for key, val in gurobi_params.items():
            if key.lower() == "logfile":
                continue
            try:
                scaled_model.setParam(key, _coerce_param_value(val))
            except gp.GurobiError as exc:
                print(
                    f"Warning: could not set Gurobi parameter "
                    f"{key}={val!r}: {exc}",
                    file=sys.stderr,
                )

        # Direct Gurobi's optimization log to our logfile.
        # optimize() will create/overwrite this file with the solve log.
        scaled_model.setParam("LogFile", logfile)

        # Solve
        scaled_model.optimize()
        opt_status = scaled_model.Status

        # Read the scaling log captured before the solve
        with open(temp_scaling_log, "r") as f:
            scaling_log_content = f.read()

        # Read Gurobi's optimization log (written by optimize())
        try:
            with open(logfile, "r") as f:
                gurobi_log_content = f.read()
        except FileNotFoundError:
            gurobi_log_content = ""

        # Compute and format quality stats if a primal solution is available
        if opt_status in (GRB.OPTIMAL, GRB.SUBOPTIMAL):
            quality_text = _format_quality_stats(scaled_model)
            print(quality_text)

        # Reassemble the integrated logfile:
        #   [scaling log] + [Gurobi solve log] + [quality stats]
        with open(logfile, "w") as f:
            f.write(scaling_log_content)
            f.write(gurobi_log_content)
            f.write(quality_text)

    finally:
        if os.path.exists(temp_scaling_log):
            os.unlink(temp_scaling_log)

    # Exit codes for non-optimal statuses
    if opt_status == GRB.INFEASIBLE:
        print("Scaled model is infeasible.", file=sys.stderr)
        sys.exit(2)
    elif opt_status == GRB.UNBOUNDED:
        print("Scaled model is unbounded.", file=sys.stderr)
        sys.exit(3)
    elif opt_status not in (GRB.OPTIMAL, GRB.SUBOPTIMAL):
        print(
            f"Solver finished with status code: {opt_status}", file=sys.stderr
        )
        sys.exit(4)
