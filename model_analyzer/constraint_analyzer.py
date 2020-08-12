from collections import defaultdict
import gurobipy as gp
import math

import model_analyzer.common as common


def get_rhs_frequencies(model, basis):
    rhs_count = defaultdict(int)
    min_exponent = float("inf")
    max_exponent = -float("inf")

    for constr in model.getConstrs():
        if constr.RHS != 0:
            exponent = int(math.floor(math.log(abs(constr.RHS), basis)))
            rhs_count[exponent] += 1
            if exponent < min_exponent: min_exponent = exponent
            if exponent > max_exponent: max_exponent = exponent

    if min_exponent == float("inf"):
        return []
    else:
        return [[exponent, rhs_count[exponent]] for exponent in range(min_exponent, max_exponent + 1)]


def get_a_frequencies(model, basis, topx):
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

                exponent = int(math.floor(math.log(abs(coeff), basis)))

                if row.getVar(i).VType == 'C':
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
                rows.append({"constr_name": constr.ConstrName,
                             "range": max_exponent - min_exponent,
                             "min_exponent": min_exponent,
                             "min_abs_coeff": min_abs_coeff,
                             "min_var": min_var.VarName,
                             "min_var_type": min_var.VType,
                             "max_exponent": max_exponent,
                             "max_abs_coeff": max_abs_coeff,
                             "max_var": max_var.VarName,
                             "max_var_type": max_var.VType
                             })

    cont_coeff_result = [[exponent, cont_coeff_count[exponent]] for exponent in
                         range(min(cont_coeff_count.keys()), max(cont_coeff_count.keys()) + 1)]
    int_coeff_result = [[exponent, int_coeff_count[exponent]] for exponent in
                        range(min(int_coeff_count.keys()), max(int_coeff_count.keys()) + 1)]
    range_result = [[exponent, range_count[exponent]] for exponent in
                    range(min(range_count.keys()), max(range_count.keys()) + 1)]

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

    data["numNZs"] = int(model.DNumNZs)

    top_ranges = []
    top_x = 10

    for i in [2, 10]:
        data["objContinuousExponentsLog" + str(i)], data["objIntegerExponentsLog" + str(i)] = \
            common.get_obj_c_frequencies(model, i)

        a_cont_coeff_count, a_int_coeff_count, a_range_count, top_ranges = get_a_frequencies(model, i, top_x)

        data["linearConstraintsContinuousExponentsLog" + str(i)] = a_cont_coeff_count
        data["linearConstraintsIntegerExponentsLog" + str(i)] = a_int_coeff_count
        data["linearConstraintsRowRangesLog" + str(i)] = a_range_count
        data["topSingleRowRangesLog" + str(i)] = top_ranges


def process_rhs(model, data):
    print("Processing RHS values...")

    for i in [2, 10]:
        data["rhsExponentsLog" + str(i)] = get_rhs_frequencies(model, i)


def process_constraints(model, data):
    print("Processing model constraints...")

    constraint_variables = set()
    redundant_constraints = []
    infeasible_constraints = []

    data["numLinConstrs"] = model.NumConstrs
    data["numQuadConstrs"] = model.NumQConstrs
    data["numSOSConstrs"] = model.NumSOS

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
                lhs_lb = max(-gp.GRB.INFINITY, lhs_lb + common.min_value(var, coefficient))

            if lhs_ub < gp.GRB.INFINITY:
                lhs_ub = min(gp.GRB.INFINITY, lhs_ub + common.max_value(var, coefficient))

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
        if (constraint.Sense == gp.GRB.LESS_EQUAL and lhs_ub <= constraint.RHS or
                constraint.Sense == gp.GRB.GREATER_EQUAL and lhs_lb >= constraint.RHS or
                constraint.Sense == gp.GRB.EQUAL and lhs_lb == constraint.RHS and lhs_ub == constraint.RHS):
            counter_redundant += 1
            redundant_constraints.append(constraint)

        # Check if constraint is infeasible
        if (constraint.Sense == gp.GRB.LESS_EQUAL and lhs_lb > constraint.RHS or
                constraint.Sense == gp.GRB.GREATER_EQUAL and lhs_ub < constraint.RHS or
                constraint.Sense == gp.GRB.EQUAL and (lhs_ub < constraint.RHS or lhs_lb > constraint.RHS)):
            counter_infeasible += 1
            infeasible_constraints.append(constraint)

        # Continue only for constraints that contain at least two variables
        if row.size() < 2: continue

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
                if sense == gp.GRB.LESS_EQUAL: sense = gp.GRB.GREATER_EQUAL
                if sense == gp.GRB.GREATER_EQUAL: sense = gp.GRB.LESS_EQUAL
            else:
                factor = 1

            # All coefficients 1?
            if abs(min_coefficient - max_coefficient) < 1e-10 and abs(min_coefficient - factor) < 1e-10:

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
            if is_row_coeff_integer and (
                    (is_row_all_positive and rhs > 0) or (is_row_all_negative and rhs < 0)) and is_rhs_integer:
                counter_bin_knapsack += 1

        # All variables integer?
        if num_intvars == row.size():

            # Knapsack constraints (all variables binary, all coefficients + RHS integer)
            if is_row_coeff_integer and (
                    (is_row_all_positive and rhs > 0) or (is_row_all_negative and rhs < 0)) and is_rhs_integer:
                counter_int_knapsack += 1

    # Loop through all quadratic constraints
    for qconstraint in model.getQConstrs():

        qcrow = model.getQCRow(qconstraint)
        row = qcrow.getLinExpr()

        # Loop quadratic terms
        for i in range(qcrow.size()):
            if qcrow.getCoeff(i) != 0:
                constraint_variables.add(qcrow.getVar1(i))
                constraint_variables.add(qcrow.getVar2(i))

        # Loop linear terms
        for i in range(row.size()):
            coeff = row.getCoeff(i)
            if coeff != 0:
                constraint_variables.add(row.getVar(i))

        if qconstraint.Sense == gp.GRB.EQUAL:
            counter_equality += 1
        else:
            counter_inequality += 1

    # TOP 5 redundant constraints
    redundant_constraints = (sorted(redundant_constraints, key=lambda constraint: model.getRow(constraint).size()))[:5]
    data["redundantConstraints"] = [constraint.ConstrName for constraint in redundant_constraints]

    # TOP 5 infeasible constraints
    infeasible_constraints = (sorted(infeasible_constraints, key=lambda constraint: model.getRow(constraint).size()))[
                             :5]
    data["infeasibleConstraints"] = [constraint.ConstrName for constraint in infeasible_constraints]

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
