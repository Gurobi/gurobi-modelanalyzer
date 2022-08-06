import gurobipy as gp
from gurobipy import GRB
from common import *
import argparse
#
#   Ill conditioning explainer.   If the model has basis statuses, it will
#   use them, computing the factorization if needed.   If no basis statuses
#   are available, solve the LP to optimality (or whatever the final status
#   is) and use that basis.
#   Computing exact Kappa can be expensive.   Set KappaExact to 0 to avoid,
#   KappaExact to 1 to force this (in addition to estimated Kappa).  Default
#   of -1 decides based on model size.
#

BASIC      =  0         # VBasis status IDs
AT_LB      = -1
SUPERBASIC = -3
SOLVELP    = 0
SOLVEQP    = 1
SOLVEMIP   = 2


def kappa_explain(model, data=None, KappaExact = -1, prmfile = None,  \
                  relaxobjtype = SOLVELP):

    if (model.IsMIP or model.IsQP or model.IsQCP):
        print("Ill Conditioning explainer only operates on LPs.")
        return None
#    import pdb; pdb.set_trace()    

    modvars = model.getVars()
    modcons = model.getConstrs()
    
    explmodel = extract_basis(model, modvars, modcons)
    explmodel.update()
    resmodel  = model.copy()
    kappastats(model, data, KappaExact)
    explmodel.write("explmodel.lp")          # debug only
#
#   Minimize the violations of B'y = 0 constraints.  Do not relax
#   the e'y == 1 constraints, and the y variables are free, so they
#   have no bounds to relax.
#
    exvars  = explmodel.getVars()
    excons  = explmodel.getConstrs()
    nbas    = len(excons) - 1
#
#   The feasRelax model involves an ill conditioned basis matrix at its
#   core, so we take steps to help the LP solve deal with that.  Start
#   with all the y variables superbasic, set NumericFocus to 3, and look at the
#   row and column ratios of the basis matrix to decide if we should
#   use geometric mean scaling.
#
#    for v in exvars:
#        v.VBasis = SUPERBASIC
    explmodel.setParam("NumericFocus", 3)
    explmodel.setParam("Method", 0)
    explmodel.setParam("FeasibilityTol", 1e-9)
    explmodel.setParam("OptimalityTol", 1e-9)
    rowrat = 1.0
    for con in excons[0:nbas]:
        rowvals = []
        rowlhs  = explmodel.getRow(con)
        for j in range(rowlhs.size()):
            rowvals.append(rowlhs.getCoeff(j))
        freqs = get_vector_frequencies(rowvals, 10)
        rat = freqs[len(freqs) - 1][0] - freqs[0][0]     # base 10 exponents
        if rat > rowrat:
            rowrat = rat
    if rowrat > 5:
        explmodel.setParam("ScaleFlag", 2)
    else:   # Row ratios don't need geometric mean scaling; check col ratios.
        colrat  = 1.0
        for var in exvars:
            colvals = []
            col = explmodel.getCol(var)
            if col.size() == 0:
                continue
            for i in range(col.size()):
                colvals.append(col.getCoeff(i))
            freqs = get_vector_frequencies(colvals, 10)  
            rat = freqs[len(freqs) - 1][0] - freqs[0][0] # base 10 exponents
            if rat > colrat:
                colrat = rat
        if colrat > 5:
            explmodel.setParam("ScaleFlag", 2)
    if prmfile != None:
        explmodel.read(prmfile)
    explmodel.update()
    explmodel.feasRelax(relaxobjtype, False, None, None, None, \
                        excons[0:nbas], [1.0]*nbas)
#
#   Solve configuration completed.  Solve the model that will give us
#   a certificate of ill conditioning.
#
    explmodel.write("explmodel_fr.mps")
    explmodel.optimize()
    status = explmodel.status
    if status in (GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED, \
                  GRB.NUMERIC, GRB.SUBOPTIMAL):
        #
        # Under perfect precision, the feasrelax model we are solving
        # should never exit infeasible or unbounded.  In that case, or if
        # the status indicates numerical issues of some sort, we should
        # take even more precautions to deal with the ill conditioned
        # submatrix in the feasrelax problem.  We do this by starting with
        # the all slack basis.   The idea is we want to find the minimal
        # relaxation with as few structural variables in the basis as
        # possible.   We don't do this by default because it was found to
        # be too time consuming.
        #
        allvars    = explmodel.getVars()
        structvars = allvars[0:nbas]
        slackvars  = allvars[nbas:]
        for struct in structvars:
            struct.VBasis = SUPERBASIC
        for slack in slackvars:
            slack.VBasis = AT_LB
        cons = explmodel.getConstrs()
        for c in cons:
            c.CBasis = BASIC
        explmodel.optimize()
        
#
#   The y variables were created before calling feasRelax.  Extract
#   them, as they are the certificate of ill conditioning.
#
    yvars     = explmodel.getVars()[0:nbas]
    rcons     = resmodel.getConstrs()
    rconsdict = {}
    count     = 0
    yvals     = []
    delcons   = []
#
#   The order in which we created variables when extracting creating 
#   the explainer problem typically will not match the order of the
#   constraints in the original model.   So we need to use a dictionary
#   to map the support of the y vector to the correct constraints
#   in the computed explanation.
#
    combinedlhs = gp.LinExpr()            # y'A; y is inf. certificate
    for c in rcons:
        rconsdict[c.ConstrName] = c
        
    for yv in yvars:
        if abs(yv.X) < 1e-13:             # TODO: make this tol relative
                                          # based on max row coeff.
            delcons.append(rconsdict[yv.VarName])  # To be filtered out.
        else:
#            print("Include constraint ", yv.VarName)     #dbg
            yvals.append(abs(yv.X))
            thiscon = rconsdict[yv.VarName]
            thiscon.ConstrName = ("(mult=" + str(yv.X) + ")") + \
                thiscon.ConstrName
            # Inf. Certificate contribution to this constraint.
            combinedlhs.add(resmodel.getRow(thiscon), yv.X)
        count += 1
                  
    resmodel.addLConstr(combinedlhs, GRB.LESS_EQUAL, GRB.INFINITY, "Combined")  
    resmodel.remove(delcons)
    resmodel.update()
#
#   Filter out the nonbasic variables of the original model that appear
#   in the constraints ih the explainer model.
#
    delvars    = []
    resvardict = {}
    resvars    = resmodel.getVars()
    for v in resvars:
        resvardict[v.VarName] = v
    for v in model.getVars():
        if v.vBasis != BASIC:
            delvars.append(resvardict[v.VarName])
    resmodel.remove(delvars)
    resmodel.update()
#
#   Print the combined value y'B to the screen, including all nonzero
#   values
#
    combinedcon = (resmodel.getConstrs())[resmodel.numConstrs - 1]
    print("Vector matrix product of certificate of ill conditioning and basis:")
    print(resmodel.getRow(combinedcon))
                           
    resmodel.setObjective(0)
    if model.ModelName == "":
        modelname = "model"
    else:
        modelname = model.ModelName
    filename = modelname + "_kappaexplain.lp"
    resmodel.write(filename)
#
#   Final info
#
    print("--------------------------------------------------------")
    print("Ill conditioning explanation written to file ", filename)
    print("Maximum absolute multiplier value: ", str(max(yvals))) 
    print("Minimum absolute multiplier value: ", str(min(yvals))) 
    print("--------------------------------------------------------")
#
#   For a given basis matrix B from the model provided, create the basis
#   model:
#   B'y = 0
#   e'y = 1         // normalization of y != 0 constraint
#   y free
#
#   This is an infeasible model for any nonsingular basis matrix B.
#
def extract_basis(model, modvars, modcons):
#
#   Does the model have a factorized basis?  If not, need to solve it
#   first.
#
    need_basis = False
    try:
        test = model.Kappa
    except AttributeError as e:
        need_basis = True     # TODO: handle statuses but no LU case.
    
    if need_basis:
        model.optimize() #TODO:check status to confirm basis available after solve.
        
    m        = model.numConstrs
    n        = model.numVars
    explmodel = gp.Model("basismodel")
#
#   Extract the basic structural variables from the model into the model
#   containing the basis matrix of interest.  Each column of the basis
#   matrix corresponds to a constraint in the explainer model.  Hence
#   each row of the basis matrix corresponds to a variable in the explainer
#   model.
#
    explvardict = {}
    for var in modvars:          # structural basic variables
        if var.VBasis != BASIC:
            continue
        col = model.getCol(var)
        varlist  = []
        coeflist = []
        for i in range(col.size()):
            coeff = col.getCoeff(i)
            con   = col.getConstr(i)
            vname = con.ConstrName
            if not vname in explvardict:
                explvardict[vname] = explmodel.addVar(lb = -float('inf'),\
                                                      name = vname)
            varlist.append(explvardict[vname])
            coeflist.append(coeff)
            
        lhs = gp.LinExpr(coeflist, varlist)
        explmodel.addConstr(lhs == 0, name=var.VarName)

    for con in modcons:          # slack basic variables
        if con.CBasis != BASIC:
            continue
        if con.Sense == '>':
            coeff = -1.0
        else:
            coeff = 1.0
        vname = con.ConstrName
        if not vname in explvardict:
            explvardict[vname] = explmodel.addVar(lb = -float('inf'),\
                                                  name = vname)
        explmodel.addConstr(coeff*explvardict[vname] == 0, \
                            name="slack_" + con.ConstrName)
#
#   B'y = 0 constraints are done.   All variables are now created, so
#   add the e'y = 1 constraint and we are done.
#
    explmodel.update()
    explmodel.addConstr(gp.quicksum(explmodel.getVars()) == 1)
    return explmodel



def kappastats(model, data, KappaExact):
    kappa      = model.Kappa
    kappaexact = -1
    if KappaExact == 1 or (model.numConstrs < 1000000 and \
                           model.numVars    < 1000000):
        kappaexact = model.KappaExact
    print("--------------------------------------------------------")
    print("Condition number stats for basis for model ", \
          model.ModelName, ":")
    print("Estimated condition number: ", kappa)
    if kappaexact == -1:
        print("Exact condition number not computed.")
    else:
        print("Exact condition number: ", kappaexact)
    print("--------------------------------------------------------")

    if data != None:
        data["Kappa"]      = kappa
        data["KappaExact"] = kappaexact
#
#   TODO: include Skeel condition number calculation
#
    
