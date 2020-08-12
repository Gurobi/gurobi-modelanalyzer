import gurobipy as gp

import model_analyzer.common as common


def process_variable_types(m, data, objective_variables, constraint_variables):
    print("Processing model variables...")

    # Count variable type
    variable_types = {gp.GRB.INTEGER: 'Gen',
                      gp.GRB.BINARY: 'Bin',
                      gp.GRB.CONTINUOUS: 'Con',
                      gp.GRB.SEMICONT: 'SemiCont',
                      gp.GRB.SEMIINT: 'SemiInt'}

    variable_count = {VType: 0 for VType in variable_types.keys()}

    variable_boundtypes = ["Free", "LBOnly", "UBOnly", "Bnd", "Fixed"]
    variable_boundcount = {BType: 0 for BType in variable_boundtypes}

    redundant_vars = []

    for var in m.getVars():

        variable_count[common.get_vtype(var)] += 1

        if var.Obj == 0:
            if m.NumPWLObjVars > 0:
                pwl = m.getPWLObj(var)
                if len(pwl) == 0:
                    column = m.getCol(var)
                    if column.size() == 0:
                        redundant_vars.append(var)

        if var.LB == -gp.GRB.INFINITY:

            if var.UB == gp.GRB.INFINITY:
                variable_boundcount["Free"] += 1
            else:
                variable_boundcount["UBOnly"] += 1

        else:

            if var.UB == gp.GRB.INFINITY:
                variable_boundcount["LBOnly"] += 1
            else:

                if var.LB == var.UB:
                    variable_boundcount["Fixed"] += 1
                else:
                    variable_boundcount["Bnd"] += 1

    for vtype, count in variable_count.items():
        data["num%sVars" % variable_types[vtype]] = count

    for btype, count in variable_boundcount.items():
        data["num%sVars" % btype] = count

    data["numVars"] = m.NumVars

    # Check for redundant variables (not in objective, not in any constraint)
    redundant_vars = list(set(m.getVars()).difference(objective_variables, constraint_variables))
    data["numRedundantVars"] = len(redundant_vars)
    data["redundantVars"] = [x.VarName for x in redundant_vars[:10]]