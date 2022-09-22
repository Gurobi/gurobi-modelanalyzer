import gurobipy as gp
from gurobipy import GRB
from common import *
import argparse
import os
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
BYROWS     = 1
BYCOLS     = 2

def kappa_explain(model, data=None, KappaExact = -1, prmfile = None,  \
                  relobjtype = SOLVELP, expltype = BYROWS, \
                  splitfreevars = True):
#
#   Help function info
#    
    '''Ill conditioning explainer.   Any basis statuses, in the model
       will be used, computing the factorization if needed.   If no statuses
       are available, solve the LP to optimality (or whatever the final 
       solve status is) and use that basis.

       Arguments:
       model      (required) The LP model whose basis will be examined.
                             Basis can be from completed solve or statuses.
                             LP will be optimized if no basis present
       data                  Do not provide.  Used only for Gurobi routines.
       KappaExact (optional) 1 = display exacti condition number in stats
                             0 = display less computationally intensive estimate
                            -1 = (default) decide based on problem size.
       prmfile    (optional) Parameter settings file for the subproblem
                             solved to derive the ill conditioning certificate.
       relobjtype (optional) Type of subproblem to create for calculation of
                             ill conditioning certificate.
                             SOLVELP (default) 
                             SOLVEQP
                             SOLVEMIP 
       expltype   (optional) Row (BYROWS) or column (BYCOLS) based computation 
                             and explanation'''
    
    if (model.IsMIP or model.IsQP or model.IsQCP):
        print("Ill Conditioning explainer only operates on LPs.")
        return None
    import pdb; pdb.set_trace()    

    modvars = model.getVars()
    modcons = model.getConstrs()
    
    explmodel, splitvardict = extract_basis(model, modvars, modcons, \
                                            expltype, splitfreevars)
    explmodel.update()
    resmodel  = None
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
    explmodel.feasRelax(relobjtype, splitfreevars, None, None, None, \
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
    yvals     = []
    yvars     = None
    if expltype == BYROWS:
        yvars     = explmodel.getVars()[0:nbas]
        resmodel  = model.copy()
        rcons     = resmodel.getConstrs()
        rconsdict = {}
        count     = 0
        delcons   = []
#
#       The order in which we created variables when extracting 
#       the explainer problem typically will not match the order of the
#       constraints in the original model.   So we need to use a dictionary
#       to map the support of the y vector to the correct constraints
#       in the computed explanation.
#

        combinedlhs = gp.LinExpr()            # y'A; y is inf. certificate
        for c in rcons:
            rconsdict[c.ConstrName] = c
        
        for yv in yvars:
            if splitfreevars:
                #
                # resmodel is in terms of original model.   Need to
                # provide y values and constraint names in terms of
                # the original model, not the explainer model in which
                # the free variables have been split.  Make use of the
                # splitvardict dictionary that connects each pair of
                # split free variables to do this.  We need to extract
                # the original model constraint name from the split
                # variable names that have the plus/minus GRB suffixes.
                #
                yval     = yv.X - splitvardict[yv.VarName].X
                ynamelen = len(yv.VarName) - len("_GRBPlus") 
                yname    = yv.VarName[0:ynamelen]   # original constr name.
            else:
                yval     = yv.X
                yname    = yv.VarName
                
            if abs(yval) < 1e-13:             # TODO: make this tol relative
                                              # based on max row coeff.
                delcons.append(rconsdict[yname])  # To be filtered out.
            else:
                print("Include constraint ", yname)     #dbg
                yvals.append(abs(yval))
                thiscon = rconsdict[yname]
                thiscon.ConstrName = ("(mult=" + str(yval) + ")") + \
                    thiscon.ConstrName
                # Inf. Certificate contribution to this constraint.
                combinedlhs.add(resmodel.getRow(thiscon), yval)
            count += 1           # do we actually need this?
                  
        resmodel.addLConstr(combinedlhs, GRB.LESS_EQUAL, GRB.INFINITY, \
                            "Combined_Row")  
        resmodel.remove(delcons)
        resmodel.update()
#
#       Filter out the nonbasic variables of the original model that appear
#       in the constraints ih the explainer model.
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
#       Print the combined value y'B to the screen, including all nonzero
#       values
#
        combinedcon = (resmodel.getConstrs())[resmodel.numConstrs - 1]
        print("Vector matrix product of certificate of ill conditioning" + \
              " and basis:")
        print(resmodel.getRow(combinedcon))
    else:              # Column based explanation
        yvars       = explmodel.getVars()[0:nbas]
        resmodel    = explmodel.copy()
        rvars       = resmodel.getVars()
        rvarsdict   = {}
        delvars     = []
        combinedcol = gp.Column()
        for v in rvars:
            rvarsdict[v.VarName] = v
        for yv in yvars:
            if splitfreevars:
                #
                # Unlike the row explanation, resmodel may include slack
                # variables that are not explicit variables of the original
                # model.   Need to provide y (i.e. column) values in terms of
                # the original model, not the explainer model in which
                # the free variables have been split.  Make use of the
                # splitvardict dictionary that connects each pair of
                # split free variables to do this.  We need to extract
                # the original model constraint name from the split
                # variable names that have the plus/minus GRB suffixes.
                #
                minusvar   = splitvardict[yv.VarName]
                yval       = yv.X - minusvar.X
                ynamelen   = len(yv.VarName) - len("_GRBPlus") 
                yname      = yv.VarName[0:ynamelen]   # original variable name
                delvars.append(rvarsdict[minusvar.VarName])
            else:
                yval     = yv.X
                yname    = yv.VarName
            if abs(yval) < 1e-13:
                # TODO: make this tol relative based on max row coeff.
                delvars.append(rvarsdict[yv.VarName])  # To be filtered out.
            else:
                #
                # Don't include relaxation variables.
                #
                print("Include variable ", yname)       #dbg
                yvals.append(abs(yval))
                thisvar = rvarsdict[yv.VarName]
                thisvar.VarName = ("(mult=" + str(yval) + ")") + yname
                thiscol   = resmodel.getCol(thisvar)
                coefflist = []
                conlist   = []
                for k in range(thiscol.size()):
                    coefflist.append(yval*thiscol.getCoeff(k))
                    conlist.append(thiscol.getConstr(k))
                combinedcol.addTerms(coefflist, conlist)

#
#       Remove the feasrelax slack variables that were introduced into the
#       explainer model.  Then remove the variables from the explainer model
#       that are not part of the explanation.
#
        explmodvars = resmodel.getVars()
        start = nbas
        if splitfreevars:
            start = 2*nbas
        end   = len(explmodvars)
        resmodel.remove(explmodvars[start:end])
        #resmodel.remove(delvars)
        for dv in delvars:
            resmodel.remove(dv)
            
        resmodel.update()
        resmodel.addVar(name="Combined_Column", column=combinedcol)
        rcons = resmodel.getConstrs()
        #
        # e'x = 1 is not in explanation
        # if using lasso approach, then feasrelax phase I constraint
        # also needs to be removed from the explanation
        #
        resmodel.remove(rcons[len(rcons) - 1])
        if splitfreevars:
            resmodel.remove(rcons[len(rcons) - 2])
        #
        # Done with column specific part of explanation.
        #
    resmodel.setObjective(0)
    if model.ModelName == "":
        modelname = "model"
    else:
        modelname = model.ModelName
    if expltype == BYROWS:
        filename = modelname + "_kappaexplain.lp"
    else:
        filename = modelname + "_kappaexplain.mps"
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
#   B'y = 0         // or By = 0 if modeltype == BYCOLS
#   e'y = 1         // normalization of y != 0 constraint
#   y free
#
#   This is an infeasible model for any nonsingular basis matrix B.
#
def extract_basis(model, modvars, modcons, modeltype=BYROWS, \
                  splitfreevars=False):
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
        #
        # No factorized basis available. If basis statuses are available,
        # use the basis specified by the statuses, irrespective of optimality.
        #
        try:
            model.write("GRBjunk.bas")
            #
            # No factorization, but complete basis statuses available.
            # Compute the explanation on the basis in the statuses.
            #
            model.setParam("IterationLimit", 0)
            model.setParam("Method", 0)
            os.remove("./GRBjunk.bas")
        except AttributeError as e:
            tmp = 0      # No statuses or factorization; solve from scratch
            
        model.optimize() #TODO:check status to confirm basis available after solve.
        
    m        = model.numConstrs
    n        = model.numVars
#
#   Extract the basic structural variables from the model into the model
#   containing the basis matrix of interest.  Each column of the basis
#   matrix corresponds to a constraint in the explainer model.  Hence
#   each row of the basis matrix corresponds to a variable in the explainer
#   model.
#
    explvardict = {}
    explmodel = gp.Model("basismodel")
    if modeltype == BYROWS:
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
    else:          # BYCOLS; assumes modeltype has been checked by caller
        explcondict = {}
        for con in modcons:
            #
            # Constraint initialization for building explainer model by col
            #
            explcondict[con.ConstrName] = \
                explmodel.addConstr(0, GRB.EQUAL, 0, name=con.ConstrName)

        #
        # Structural basic variables.  The column from the original model
        # needs to reference the constraint in the explainer model, not
        # the original model.
        #
        for var in modvars:          
            if var.VBasis != BASIC:
                continue
            col = model.getCol(var)
            for k in range(col.size()):
                colcon = col.getConstr(k)
                colcon = explcondict[colcon.ConstrName]
            explmodel.addVar(obj=0.0, lb = -float('inf'), name=var.VarName, \
                             column=col)
        #
        #   Structural variables done.  Now do the basic slack variables
        #
        for con in modcons:          # slack basic variables
            if con.CBasis != BASIC:
                continue
            if con.Sense == '>':
                coeff = -1.0
            else:
                coeff = 1.0
            vname = "slack_" + con.ConstrName
            col = gp.Column([coeff], [explcondict[con.ConstrName]])
            explmodel.addVar(obj=0.0, lb = -float('inf'), name=vname, \
                             column=col)
#
#   B'y = 0 or By = 0 constraints are done.   All variables are now created, so
#   add the e'y = 1 constraint and we are done.  This is the same,
#   regardless of whether the explanation model is row or column based.
#
    explmodel.update()
    splitvardict = None
    if splitfreevars:
        plusvars     = explmodel.getVars()
        splitvardict = splitmirroredvars(explmodel)
        minusvars    = []
        for v in plusvars:
            minusvar = splitvardict[v.VarName]
            minusvars.append(minusvar)
#            explmodel.addSOS(GRB.SOS_TYPE1, [v, minusvar]) 
                             
        explmodel.addConstr(gp.quicksum(plusvars) - gp.quicksum(minusvars) == 1)
        explmodel.setObjective(gp.quicksum(plusvars) + gp.quicksum(minusvars))
    else:
        explmodel.addConstr(gp.quicksum(explmodel.getVars()) == 1)        

    return explmodel,splitvardict



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
#   Split free variables into difference of two nonnegative variables.
#   Also works on other variables whose lower bound is the opposite of
#   its upper bound.  varstosplit identifies the candidate list of variables
#   to split.  If not set or None, will treat all variables as candidates.
#   Only those that meet the "mirrored variable" criteria will
#   be split.   Most common use case will be to split free variables.
#
#   Alters the model, so use a copy if you want to preserve the original
#   model.
#
#   Tests:  1) Test on a model with only QCs.
#           2) Test on a model with linear obj terms but no linear constraints.
#
def splitmirroredvars(model, varstosplit=None):
    model.update()
    if varstosplit == None:
        varstosplit = model.getVars()

    newvarlist = []           # This var and dict is used to handle
    newvardict = {}           # quadratic objective and constraints
    for v in varstosplit:
        if v.UB != -v.LB:
            continue
        #
        #  We have a candidate variable to split.
        #
        bnd        = v.UB
        varname    = v.VarName 
        chgvarname = v.VarName  + "_GRBPlus"
        v.VarName  = chgvarname
        col        = model.getCol(v)
        splitcol   = gp.Column()
        coeflist   = []
        conlist    = []
        for k in range(col.size()):
            coeflist.append(-col.getCoeff(k))
            conlist.append(col.getConstr(k))
        splitcol.addTerms(coeflist, conlist)
        v.LB = 0
        newvar = model.addVar(lb=0.0, ub=bnd, obj=-v.obj, vtype=v.Vtype, \
                              name=varname + "_GRBMinus", column=splitcol)
        newvarlist.append(newvar)
        newvardict[chgvarname] = newvar

    #
    # Linear constraints completed; now update any QCs that contain
    # mirrored variables that need to be split.
    #
    for qcon in model.getQConstrs():
        quadexpr = model.getQCRow(qcon)
        splitquadexpr(quadexpr, newvardict)
        #
        # Quadratic expression is updated; need to create a new
        # QC with the update quadexpr and delete the old one.
        #
        model.addQConstr(quadexpr, qcon.QCSense, qcon.QCRHS, qcon.QCName)
        model.remove(qcon)
                    
    #
    # Quadratic constraints completed; now the quadratic objective if it 
    # contains any mirrored variables that need to be split.
    #
    objexpr = model.getObjective()
    if isinstance(objexpr, gp.QuadExpr):
        splitquadexpr(objexpr, newvardict)
        model.setObjective(objexpr)
    model.update()
    return newvardict
#
#   Takes a quadratic expression and replaces mirrored variables with
#   the appropriate difference of two nonnegative variables.  All
#   mirrored variables are provided in newvardict, which maps the
#   original mirrored variable to its added variable.
    
def splitquadexpr(quadexpr, newvardict):
    for k in range(quadexpr.size()):
        xi = quadexpr.getVar1(k)
        xj = quadexpr.getVar2(k)
        qcoef = quadexpr.getCoeff(k)
        if xi.varName in newvardict:      # xi is split into xi(+) - xi(-)
            # xi(-) * xj term
            quadexpr.addTerms(-qcoef, xj, newvardict[xi.VarName])
            if xj.varName in newvardict: # both xi and xj split
                # xi * xj(-) term
                quadexpr.addTerms(-qcoef, xi, newvardict[xj.VarName])
                # xi(-) * xj(-) term
                quadexpr.addTerms(qcoef, newvardict[xi.VarName], \
                                      newvardict[xj.VarName])
        else:                             
            if xj.varName in newvardict: # xj is split, but xi is not.
                # xi * xj(-) term
                quadexpr.addTerms(-qcoef, xi, newvardict[xj.VarName])




                
#
#   TODO: include Skeel condition number calculation
#
    