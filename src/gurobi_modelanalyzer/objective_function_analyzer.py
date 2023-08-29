from collections import defaultdict
import gurobipy as gp

import gurobi_modelanalyzer.common as common


def process_objective(m, data):
    print("Processing objective function...")

    obj = m.getObjective()
    obj_variables = set()
    count_obj_coeffs_decimal_digits = defaultdict(int)
    count_positive_coeffs = 0
    count_linear_vars = 0
    count_pwl_vars = 0
    count_pwl_convex_vars = 0
    count_pwl_nonconvex_vars = 0
    count_quad_terms = 0

    # Quadratic terms
    if type(obj) is gp.QuadExpr:
        count_quad_terms = obj.size()
        for i in range(obj.size()):
            coeff = obj.getCoeff(i)
            if coeff > 0:
                count_positive_coeffs += 1
            count_obj_coeffs_decimal_digits[common.count_decimal_digits(coeff)] += 1
            obj_variables.add(obj.getVar1(i))
            obj_variables.add(obj.getVar2(i))
            # TODO: Include bounds for quadratic terms in obj_lb and obj_ub

        obj = obj.getLinExpr()

    count_quad_vars = len(obj_variables)
    obj_lb = obj.getConstant()
    obj_ub = obj.getConstant()

    # Linear / PWL terms
    for var in m.getVars():
        if var.Obj == 0:
            pwl = m.getPWLObj(var)
            if len(pwl) > 0:
                count_pwl_vars += 1
                obj_variables.add(var)
                if var.PWLObjCvx == 0:
                    count_pwl_nonconvex_vars += 1
                else:
                    count_pwl_convex_vars += 1
        else:
            count_linear_vars += 1
            obj_variables.add(var)
            if var.Obj > 0:
                count_positive_coeffs += 1

            if obj_lb > -gp.GRB.INFINITY:
                obj_lb = max(-gp.GRB.INFINITY, obj_lb + common.min_value(var, var.Obj))

            if obj_ub < gp.GRB.INFINITY:
                obj_ub = min(gp.GRB.INFINITY, obj_ub + common.max_value(var, var.Obj))

            count_obj_coeffs_decimal_digits[common.count_decimal_digits(var.Obj)] += 1

    obj_coeffs_decimal_digits_result = [
        [length, count_obj_coeffs_decimal_digits[length]]
        for length in range(
            min(count_obj_coeffs_decimal_digits.keys()),
            max(count_obj_coeffs_decimal_digits.keys()) + 1,
        )
    ]

    data["objLB"] = obj_lb
    data["objUB"] = obj_ub

    data["objDecimalDigits"] = obj_coeffs_decimal_digits_result
    data["numObjTotalVars"] = len(obj_variables)
    data["numObjVars"] = count_linear_vars
    data["numPWLObj"] = count_pwl_vars
    data["numPWLObjCVars"] = count_pwl_convex_vars
    data["numPWLObjNCVars"] = count_pwl_nonconvex_vars
    data["numQObjTerms"] = count_quad_terms
    data["numQObjVars"] = count_quad_vars

    data["numObjIntVars"] = len(
        [var for var in obj_variables if common.get_vtype(var) != gp.GRB.CONTINUOUS]
    )
    data["numObjCoeffs"] = sum(count_obj_coeffs_decimal_digits.values())
    data["numObjIntCoeffs"] = count_obj_coeffs_decimal_digits[0]
    data["numObjPosCoeffs"] = count_positive_coeffs

    if m.ModelSense == 1:
        data["isMin"] = True
    else:
        data["isMin"] = False

    assert count_pwl_vars == count_pwl_convex_vars + count_pwl_nonconvex_vars
    return obj_variables
