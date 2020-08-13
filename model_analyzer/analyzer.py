import contextlib
import datetime
import gurobipy as gp
from io import StringIO
import os
import sys
from timeit import default_timer

import model_analyzer.common as common
import model_analyzer.constraint_analyzer as constraint_analyzer
import model_analyzer.file_analyzer as file_analyzer
import model_analyzer.objective_function_analyzer as objective_function_analyzer
import model_analyzer.variable_analyzer as variable_analyzer


def run(model_file: str, analyzers):
    # Uncompress model file if compressed
    fileType, compression = file_analyzer.get_file_type(model_file)

    if fileType is None:
        print(f"Error: {model_file} has unsupported file type")
        return 1

    # Initialize data dictionary
    data = {
        "fileName": os.path.basename(model_file),
        "fileType": fileType.upper(),
        "reportTimestamp": common.utc_date_formatted(datetime.datetime.utcnow()),
        "gurobiVersion": "%d.%d.%d" % gp.gurobi.version(),
        "fileSize": os.stat(model_file).st_size,
        "fileSizeUncompressed": 0,
        "readerSuccess": False,
        "readerWarnings": False,
        "readerError": False
    }

    # Read model into Gurobi
    try:
        print(f"Reading model file '{model_file}' ...")
        with elapsed_timer() as elapsed:
            with capture_output() as out:
                m = gp.read(model_file)
    except:
        data["readerError"] = True

    finally:
        data["readerOutput"] = [line.rstrip('\n') for line in out if len(line) > 0]
        data["readerTime"] = elapsed()

    if data["readerError"]:
        print("ERROR reading model file")
        return 1

    if compression is not None:

        data["isCompressed"] = True

        try:
            file_analyzer.process_compressed_model(compression, model_file, data)

        except Exception as e:
            print("ERROR while trying do decompress model file")
            print(e)
            return 1
    else:

        data["isCompressed"] = False
        file_analyzer.process_uncompressed_model(model_file, data)

    # Check for warnings
    if any("Warning" in line for line in data["readerOutput"]):
        data["readerWarnings"] = True
    else:
        data["readerSuccess"] = True

    process_model_type(m, data)
    if "objective" in analyzers:
        objective_variables = objective_function_analyzer.process_objective(m, data)

    if "constraint" in analyzers:
        constraint_variables = constraint_analyzer.process_constraints(m, data)

    if "variable" in analyzers:
        variable_analyzer.process_variable_types(m, data, objective_variables, constraint_variables)

    if "bounds" in analyzers:
        constraint_analyzer.process_bounds(m, data)

    if "coefficients" in analyzers:
        constraint_analyzer.process_coefficients(m, data)

    if "rhs" in analyzers:
        constraint_analyzer.process_rhs(m, data)

    # Write matrix plot image and put final dimensions into data
    #write_matrix_plot(m, 1920, 1080, args.outputprefix + ".png", data)

    return data


def process_model_type(m, data):
    data["isMIP"] = bool(m.IsMIP)
    data["isQP"] = bool(m.IsQP)
    data["isQCP"] = bool(m.IsQCP)

    part1_short = "MI" if data["isMIP"] == 1 else ""
    part1_long = "Mixed-Integer " if data["isMIP"] == 1 else ""

    part2_short = "Q" if data["isQP"] == 1 or data["isQCP"] == 1 else ""
    part2_long = "Quadratic " if data["isQP"] and not data["isQCP"] else ""

    part3_short = "C" if data["isQCP"] else ""
    part3_long = "Quadratically Constrained " if data["isQCP"] else ""

    part4_short = "LP" if part1_short == "" and part2_short == "" and part3_short == "" else "P"
    part4_long = "Linear Program" if part1_long == "" and part2_long == "" and part3_long == "" else "Program"

    data["modelTypeShort"] = part1_short + part2_short + part3_short + part4_short
    data["modelTypeLong"] = part1_long + part2_long + part3_long + part4_long
    print("Model Type: %s" % data["modelTypeShort"])


@contextlib.contextmanager
def capture_output():
    oldout,olderr = sys.stdout, sys.stderr
    try:
        out=[StringIO(), StringIO()]
        sys.stdout,sys.stderr = out
        yield out
    finally:
        sys.stdout,sys.stderr = oldout, olderr
        out[0] = out[0].getvalue()
        out[1] = out[1].getvalue()


@contextlib.contextmanager
def elapsed_timer():
    start = default_timer()
    elapser = lambda: default_timer() - start
    yield lambda: elapser()
    end = default_timer()
    elapser = lambda: end-start