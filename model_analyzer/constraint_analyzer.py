from collections import defaultdict
import gurobipy as gp
import math

from model_analyzer import common


def get_rhs_frequencies(model, base=10):
    rhs_count = defaultdict(int)
    rhs_vals = []
    result = []

    for constr in model.getConstrs():
        if constr.RHS != 0:
            rhs_vals.append(constr.RHS)

    result = common.get_vector_frequencies(rhs_vals, base)

    return result


def get_a_frequencies(model, base, topx):
    int_coeff_count = defaultdict(int)
    cont_coeff_count = defaultdict(int)
    range_count = defaultdict(int)
    rows = []

    for constr in model.getConstrs():
        row = model.getRow(constr)

        min_exponent = float("inf")
        min_abs_coeff = float("inf")
        max_exponent = -float("inf")
        max_abs_coeff = -float("inf")

        # Loop through all terms, count single coefficient exponents and total range row
        for i in range(row.size()):
            coeff = row.getCoeff(i)
            if coeff != 0:
                exponent = int(math.floor(math.log(abs(coeff), base)))

                if row.getVar(i).VType == "C":
                    cont_coeff_count[exponent] += 1
                else:
                    int_coeff_count[exponent] += 1

                if exponent < min_exponent:
                    min_exponent = exponent
                    min_var = row.getVar(i)
                    min_abs_coeff = abs(coeff)
                if exponent > max_exponent:
                    max_exponent = exponent
                    max_var = row.getVar(i)
                    max_abs_coeff = abs(coeff)

        # Record range for non-empty constraints
        if min_exponent < float("inf"):
            range_count[max_exponent - min_exponent] += 1

            if topx > 0:
                rows.append(
                    {
                        "constr_name": constr.ConstrName,
                        "range": max_exponent - min_exponent,
                        "min_exponent": min_exponent,
                        "min_abs_coeff": min_abs_coeff,
                        "min_var": min_var.VarName,
                        "min_var_type": min_var.VType,
                        "max_exponent": max_exponent,
                        "max_abs_coeff": max_abs_coeff,
                        "max_var": max_var.VarName,
                        "max_var_type": max_var.VType,
                    }
                )

    cont_coeff_result = [
        [exponent, cont_coeff_count[exponent]]
        for exponent in range(
            min(cont_coeff_count.keys()), max(cont_coeff_count.keys()) + 1
        )
    ]
    int_coeff_result = [
        [exponent, int_coeff_count[exponent]]
        for exponent in range(
            min(int_coeff_count.keys()), max(int_coeff_count.keys()) + 1
        )
    ]
    range_result = [
        [exponent, range_count[exponent]]
        for exponent in range(min(range_count.keys()), max(range_count.keys()) + 1)
    ]

    # Get the TOP x results
    if topx > 0:
        topx_results = sorted(rows, key=lambda k: k["range"], reverse=True)[:topx]
    else:
        topx_results = []

    return cont_coeff_result, int_coeff_result, range_result, topx_results


def process_bounds(model, data):
    print("Processing bounds...")

    for i in [2, 10]:
        lb, ub = common.get_variable_bounds(model, data, i)
        data["varLowerBoundsLog" + str(i)] = lb
        data["varUpperBoundsLog" + str(i)] = ub


def process_coefficients(model, data):
    print("Processing coefficients...")

    top_ranges = []
    top_x = 10

    for i in [2, 10]:
        (
            data["objContinuousExponentsLog" + str(i)],
            data["objIntegerExponentsLog" + str(i)],
        ) = common.get_obj_c_frequencies(model, i)

        (
            a_cont_coeff_count,
            a_int_coeff_count,
            a_range_count,
            top_ranges,
        ) = get_a_frequencies(model, i, top_x)

        data["linearConstraintsContinuousExponentsLog" + str(i)] = a_cont_coeff_count
        data["linearConstraintsIntegerExponentsLog" + str(i)] = a_int_coeff_count
        data["linearConstraintsRowRangesLog" + str(i)] = a_range_count
        data["topSingleRowRangesLog" + str(i)] = top_ranges


def process_rhs(model, data):
    print("Processing RHS values...")

    for i in [2, 10]:
        data["rhsExponentsLog" + str(i)] = get_rhs_frequencies(model, i)


def process_linear_constraints(model, data):
    print("Processing model linear constraints...")

    constraint_variables = set()
    redundant_constraints = []
    infeasible_constraints = []

    # Equality constraint
    counter_equality = 0

    # Inequality constraint
    counter_inequality = 0

    # Empty constraint
    counter_empty = 0

    # Bound constraint
    counter_bound = 0

    # Integer and continuous variables
    counter_mixed = 0

    # Only (>= 2) continuous variables
    counter_con_only = 0

    # Only (>= 2) binary variables
    counter_bin_only = 0

    # Only (>= 2) integer variables (including binaries)
    counter_int_only = 0

    # sum xi >= 1 (xi bin)
    counter_set_covering = 0

    # sum xi = 1 (xi bin)
    counter_set_partitioning = 0

    # sum xi <= 1 (xi bin)
    counter_set_packing = 0

    # sum xi <= c (xi bin, c integer != 1)
    counter_max_binsum = 0

    # sum xi = c (xi bin, c integer != 1)
    counter_binsum = 0

    # sum xi >= c (xi bin, c integer != 1)
    counter_min_binsum = 0

    # sum c_i x_i <= d (xi bin, ci, d positive integer)
    counter_bin_knapsack = 0

    # sum c_i x_i <= d (xi integer, ci, d positive integer)
    counter_int_knapsack = 0

    # Redundant constraints (i.e. constraint is always satisfied)
    counter_redundant = 0

    # Infeasible constraints (i.e. constraint can never be satisfied)
    counter_infeasible = 0

    # Loop through all linear constraints
    for constraint in model.getConstrs():
        if constraint.Sense == gp.GRB.EQUAL:
            counter_equality += 1
        else:
            counter_inequality += 1

        row = model.getRow(constraint)

        rhs = constraint.RHS - row.getConstant()
        is_rhs_integer = abs(rhs - round(rhs)) < 1e-10

        # Empty constraints
        if row.size() == 0:
            counter_empty += 1

        # Empty or bound constraint
        elif row.size() == 1:
            counter_bound += 1

        min_coefficient = gp.GRB.INFINITY
        max_coefficient = -gp.GRB.INFINITY

        num_binvars = 0
        num_intvars = 0

        is_row_coeff_integer = True
        is_row_all_positive = True
        is_row_all_negative = True

        lhs_lb = row.getConstant()
        lhs_ub = row.getConstant()

        for i in range(row.size()):
            var = row.getVar(i)
            coefficient = row.getCoeff(i)

            if lhs_lb > -gp.GRB.INFINITY:
                lhs_lb = max(
                    -gp.GRB.INFINITY, lhs_lb + common.min_value(var, coefficient)
                )

            if lhs_ub < gp.GRB.INFINITY:
                lhs_ub = min(
                    gp.GRB.INFINITY, lhs_ub + common.max_value(var, coefficient)
                )

            if coefficient != 0:
                constraint_variables.add(var)

            min_coefficient = min(min_coefficient, coefficient)
            max_coefficient = max(max_coefficient, coefficient)

            if (coefficient - round(coefficient)) > 1e-10:
                is_row_coeff_integer = False

            if coefficient < 0:
                is_row_all_positive = False

            if coefficient > 0:
                is_row_all_negative = False

            if common.get_vtype(var) == gp.GRB.BINARY:
                num_binvars += 1
            elif common.get_vtype(var) == gp.GRB.INTEGER:
                num_intvars += 1

        # Check if constraint is redundant
        if (
            constraint.Sense == gp.GRB.LESS_EQUAL
            and lhs_ub <= constraint.RHS
            or constraint.Sense == gp.GRB.GREATER_EQUAL
            and lhs_lb >= constraint.RHS
            or constraint.Sense == gp.GRB.EQUAL
            and lhs_lb == constraint.RHS
            and lhs_ub == constraint.RHS
        ):
            counter_redundant += 1
            redundant_constraints.append(constraint)

        # Check if constraint is infeasible
        if (
            constraint.Sense == gp.GRB.LESS_EQUAL
            and lhs_lb > constraint.RHS
            or constraint.Sense == gp.GRB.GREATER_EQUAL
            and lhs_ub < constraint.RHS
            or constraint.Sense == gp.GRB.EQUAL
            and (lhs_ub < constraint.RHS or lhs_lb > constraint.RHS)
        ):
            counter_infeasible += 1
            infeasible_constraints.append(constraint)

        # Continue only for constraints that contain at least two variables
        if row.size() < 2:
            continue

        # Int-Only / Cont-Only / Mixed
        if num_binvars + num_intvars == 0:
            counter_con_only += 1
        elif num_binvars + num_intvars == row.size():
            counter_int_only += 1
        else:
            counter_mixed += 1

        # All variables binary?
        if num_binvars == row.size():
            counter_bin_only += 1

            sense = constraint.Sense

            if rhs < 0:
                factor = -1
                if sense == gp.GRB.LESS_EQUAL:
                    sense = gp.GRB.GREATER_EQUAL
                if sense == gp.GRB.GREATER_EQUAL:
                    sense = gp.GRB.LESS_EQUAL
            else:
                factor = 1

            # All coefficients 1?
            if (
                abs(min_coefficient - max_coefficient) < 1e-10
                and abs(min_coefficient - factor) < 1e-10
            ):
                # RHS = 1?
                if abs(rhs - factor) < 1e-10:
                    if constraint.Sense == gp.GRB.EQUAL:
                        counter_set_partitioning += 1
                    elif constraint.Sense == gp.GRB.LESS_EQUAL:
                        counter_set_packing += 1
                    else:
                        counter_set_covering += 1
                    continue

                # RHS integer?
                elif is_rhs_integer:
                    if constraint.Sense == gp.GRB.EQUAL:
                        counter_binsum += 1
                    elif constraint.Sense == gp.GRB.LESS_EQUAL:
                        counter_max_binsum += 1
                    else:
                        counter_min_binsum += 1
                    continue

            # Binary knapsack constraints (all variables binary, all coefficients + RHS integer)
            if (
                is_row_coeff_integer
                and (
                    (is_row_all_positive and rhs > 0)
                    or (is_row_all_negative and rhs < 0)
                )
                and is_rhs_integer
            ):
                counter_bin_knapsack += 1

        # All variables integer?
        if num_intvars == row.size():
            # Knapsack constraints (all variables binary, all coefficients + RHS integer)
            if (
                is_row_coeff_integer
                and (
                    (is_row_all_positive and rhs > 0)
                    or (is_row_all_negative and rhs < 0)
                )
                and is_rhs_integer
            ):
                counter_int_knapsack += 1

    # TOP 5 redundant constraints
    redundant_constraints = (
        sorted(
            redundant_constraints,
            key=lambda constraint: model.getRow(constraint).size(),
        )
    )[:5]
    data["redundantConstraints"] = [
        constraint.ConstrName for constraint in redundant_constraints
    ]

    # TOP 5 infeasible constraints
    infeasible_constraints = (
        sorted(
            infeasible_constraints,
            key=lambda constraint: model.getRow(constraint).size(),
        )
    )[:5]
    data["infeasibleConstraints"] = [
        constraint.ConstrName for constraint in infeasible_constraints
    ]

    data["numEqConstrs"] = counter_equality
    data["numIneqConstrs"] = counter_inequality

    data["numEmptyConstrs"] = counter_empty
    data["numBndConstrs"] = counter_bound
    data["numRedundantConstrs"] = counter_redundant
    data["numInfeasibleConstrs"] = counter_infeasible

    data["numContConstrs"] = counter_con_only
    data["numIntConstrs"] = counter_int_only
    data["numBinConstrs"] = counter_bin_only
    data["numMixedConstrs"] = counter_mixed

    data["numSetCovConstrs"] = counter_set_covering
    data["numSetPartConstrs"] = counter_set_partitioning
    data["numSetPackConstrs"] = counter_set_packing

    data["numMaxBinSumConstrs"] = counter_max_binsum
    data["numMinBinSumConstrs"] = counter_min_binsum
    data["numBinSumConstrs"] = counter_binsum

    data["numBinKnapsackConstrs"] = counter_bin_knapsack
    data["numIntKnapsackConstrs"] = counter_int_knapsack

    return constraint_variables


def process_quadratic_constraints(model, data):
    #
    #   TODO: identify infeasible and redundant QCs.
    #   Not sure worth doing for QCs, but what about semi cont. and
    #   semi int variables?
    #

    print("Processing model quadratic constraints...")

    quadpart_variables = set()
    linpart_variables = set()

    redundant_constraints = []
    infeasible_constraints = []

    # Equality constraint
    counter_equality = 0

    # Inequality constraint
    counter_inequality = 0

    # Empty constraint
    counter_empty = 0

    # Bound constraint
    counter_bound = 0

    # Integer and continuous variables
    if model.NumIntVars + model.NumBinVars >= 1:
        counter_mixed = 0
    else:
        counter_mixed = -1

    # Only (>= 2) continuous variables

    if model.NumVars - model.NumIntVars - model.NumBinVars >= 2:
        counter_cont_only = 0
    else:
        counter_cont_only = -1  # Can't have bilinear terms with 2 cont. vars

    # Only (>= 2) binary variables
    if model.NumBinVars >= 2:
        counter_bin_only = 0
    else:
        counter_bin_only = -1  # Can't have any binary bilinear terms

    # Only (>= 2) integer variables (including binaries)
    if model.NumIntVars >= 1:
        counter_int_only = 0
    else:
        counter_int_only = -1  # Can't have any integer bilinear terms

    # No linear terms in the QC (bot not totally empty)

    counter_quad_only = 0

    # Empty quadratic part; should have made it a linear constraint

    counter_lin_only = 0

    # At least one binary in every quadratic term of the QC

    if model.NumBinVars >= 1:
        counter_common_bin = 0
    else:
        counter_common_bin = -1

    # Anything else regarding QC term composition

    counter_other = 0

    # Loop through all quadratic constraints

    BIN_ONLY = 1
    INT_ONLY = 2
    CONT_ONLY = 4
    MIXED = 8
    for qconstraint in model.getQConstrs():
        qcrow = model.getQCRow(qconstraint)
        row = qcrow.getLinExpr()
        qlen = qcrow.size()
        linlen = row.size()
        bittype = 0
        qclhs_lb = row.getConstant()
        qclhs_ub = row.getConstant()

        # Loop quadratic terms

        if qlen == 1:  #  Check for QC singleton that is a simple bound
            if qcrow.getCoeff(0) != 0:
                v1 = qcrow.getVar1(0)
                v2 = qcrow.getVar2(0)
                if v1.VarName == v2.VarName and linlen == 0:
                    counter_bound += 1
                if v1.vType == gp.GRB.BINARY:
                    counter_common_bin += 1
                    if v2.vType == gp.GRB.BINARY:
                        bittype != BIN_ONLY
                    else:
                        bittype != MIXED
                elif v2.vType == gp.GRB.BINARY:
                    counter_common_bin += 1
                    bittype != MIXED
                else:  # QC singleton of two continuous variables
                    bittype != CONT_ONLY
            elif linlen == 0:
                counter_empty += 1
            else:  #  Empty quad part with nonempty linear part
                counter_lin_only += 1
        else:  # Multiple QC terms; count all pairwise combos
            binvardict = defaultdict(int)  # for common binary factors
            for j in range(qlen):
                if qcrow.getCoeff(j) != 0:
                    v1 = qcrow.getVar1(j)
                    v2 = qcrow.getVar2(j)
                    quadpart_variables.add(v1)
                    quadpart_variables.add(v2)
                    if counter_cont_only >= 0:
                        if v1.Vtype == gp.GRB.CONTINUOUS:
                            if v2.VType == gp.GRB.CONTINUOUS:
                                bittype |= CONT_ONLY
                            else:
                                bittype |= MIXED
                                if v2.Vtype == gp.GRB.BINARY:
                                    binvardict[v2.VarName] += 1
                            continue
                        elif v2.Vtype == gp.GRB.CONTINUOUS:
                            bittype |= MIXED
                            if v2.Vtype == gp.GRB.BINARY:
                                binvardict[v2.VarName] += 1
                            continue
                    if counter_bin_only >= 0:
                        if v1.Vtype == gp.GRB.BINARY and v2.VType == gp.GRB.BINARY:
                            bittype |= BIN_ONLY
                            binvardict[v1.VarName] += 1
                            if v1.VarName != v2.VarName:
                                binvardict[v1.VarName] += 1
                            continue
                    if counter_int_only >= 0:
                        if v1.Vtype == gp.GRB.INTEGER and v2.VType == gp.GRB.INTEGER:
                            bittype |= INT_ONLY
                            continue
            if len(binvardict) == qlen:
                counter_common_bin += 1

        # Finished processing this QC; update the counts of the different
        # types of bilinear terms

        if bittype == BIN_ONLY:
            counter_bin_only += 1
        elif bittype == INT_ONLY:
            counter_int_only += 1
        elif bittype == CONT_ONLY:
            counter_cont_only += 1
        elif bittype == MIXED:
            counter_mixed += 1
        else:  # Anything that doesn't fit one of the above QC types
            counter_other += 1

        # Loop linear terms.   Record contribution of linear part of QC
        # to the QC inf and sup (any constant term already recorded at
        # initialization of inf and sup).

        for i in range(linlen):
            coeff = row.getCoeff(i)
            if coeff != 0:
                var = row.getVar(i)
                linpart_variables.add(var)
                if qclhs_lb > -gp.GRB.INFINITY:
                    t = common.min_value(var, coeff)
                    if t != -gp.GRB.INFINITY:
                        qclhs_lb += t
                    else:
                        qclhs_lb = -gp.GRB.INFINITY

                if qclhs_ub < gp.GRB.INFINITY:
                    t = common.max_value(var, coeff)
                    if t != gp.GRB.INFINITY:
                        qclhs_ub += t
                    else:
                        qclhs_ub = gp.GRB.INFINITY

        if qconstraint.QCSense == gp.GRB.EQUAL:
            counter_equality += 1
        else:
            counter_inequality += 1

        if qlen == 0:
            if linlen == 0:
                counter_empty += 1
            else:
                counter_lin_only += 1
        else:
            if linlen == 0:
                counter_quad_only += 1
    #
    #   Final results.  Counters that are at -1 were flagged as
    #   not needing counting due to the numbers of different types of
    #   variables in the model (e.g. you can't have constraints involving
    #   bilinear terms of continuous variables when the model is all binary,
    #   so flag that and don't bother counting in the inner loop across the QC
    #
    data["numQCEqConstraints"] = counter_equality
    data["numQCInEqConstraints"] = counter_inequality
    data["numQCEmptyConstrs"] = counter_empty
    data["numQCBndConstrs"] = counter_bound
    if counter_cont_only == -1:
        counter_cont_only = 0
    data["numQCContConstrs"] = counter_cont_only
    data["numQCIntConstrs"] = counter_int_only
    if counter_bin_only == -1:
        counter_bin_only = 0
    data["numQCBinConstrs"] = counter_bin_only
    if counter_mixed == -1:
        counter_mixed = 0
    data["numQCMixedConstrs"] = counter_mixed
    data["numQCOtherConstrs"] = counter_other
    data["numQCQuadOnlyConstrs"] = counter_quad_only
    data["numQCLinOnlyConstrs"] = counter_lin_only
    data["numQCEmptyConstrs"] = counter_empty
    if counter_common_bin == -1:
        counter_common_bin = 0
    data["numQCBinFactorConstrs"] = counter_common_bin

    return quadpart_variables, linpart_variables
