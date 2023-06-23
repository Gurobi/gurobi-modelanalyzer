import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys
import os
import math
import time

from model_analyzer import common


#
#   Ill conditioning explainer.   If the model has basis statuses, it will
#   use them, computing the factorization if needed.   If no basis statuses
#   are available, solve the LP to optimality (or whatever the final status
#   is) and use that basis.
#   Computing exact Kappa can be expensive.   Set KappaExact to 0 to avoid,
#   KappaExact to 1 to force this (in addition to estimated Kappa).  Default
#   of -1 decides based on model size.
#
OFF         = 0
MODERATE    = 1
VERBOSE     = 2
_debug      = OFF            # Change to MODERATE or VERBOSE as needed
_debugger   = OFF

SOLVELP     = 0              # relobjtype choices
SOLVEQP     = 1
SOLVEMIP    = 2
BYROWS      = 1              # expltype choices
BYCOLS      = 2
DEFAULT     = 0              # method choices.  Default = no regularization.
ANGLES      = 1              
LASSO       = 2              # One norm regularization.
RLS         = 3              # Two norm regularization. 
CON         = 0
VAR         = 1
COEFF       = 2
DEFSMALLTOL = 1e-13      
COMBINEDROW = "GRB_Combined_Row"
COMBINEDCOL = "GRB_Combined_Column"

def kappa_explain(model, data=None, KappaExact=-1, prmfile=None,  \
                  relobjtype="LP", expltype="ROWS", method="DEFAULT", \
                  smalltol=DEFSMALLTOL, submatrix=False):
#
#   Help function info  
#    
    '''Ill conditioning explainer.   Any basis statuses, in the model
       will be used, computing the factorization if needed.   If no statuses
       are available, solve the LP to optimality (or whatever the final 
       solve status is) and use that basis.

       Arguments:
       model      (required) The LP model whose basis will be examined.
                             Basis can be from a completed solve or statuses.
                             LP will be optimized if no basis present
       data                  Do not provide.  Used only for Gurobi routines.
       KappaExact (optional) 1 = display exact condition number in stats
                             0 = display less computationally intensive estimate
                            -1 = (default) decide based on problem size.
       prmfile    (optional) Parameter settings file for the subproblem
                             solved to derive the ill conditioning certificate.
       relobjtype (optional) Type of subproblem to create for calculation of
                             ill conditioning certificate.
                             "LP" (default) 
                             "QP"
       expltype   (optional) Row ("ROWS") or column ("COLS") based computation 
                             and explanation. Default is by row.  One type of
                             explanation is often much smaller than the other.
       method     (optional) Alternate computation methods (may help generate
                             smaller explanation).
                             "DEFAULT" = Basic subproblem (no regularization). 
                             "ANGLES"  = Simpler calculation based on inner 
                                         products.  Potentially faster but only 
                                         finds explanations of 2 rows or 
                                         columns.
                             "LASSO"   = One norm regularization.
                             "RLS"     = Two norm regularization.
       smalltol   (optional) Tolerance below which certificate of ill 
                             conditioning values are treated as zero.   If
                             left at default of 1e-13, row or column norm and
                             machine precision will be incorporated.
       submatrix  (optional) Whether to postprocess the explanation down to 
                             a smaller submatrix.  Default is False.
       Returns:              For all method settings except ANGLES, returns
                             the model that consists of the explanation.  
                             If the method is ANGLES, returns a list of tuples 
                             of the almost parallel rows and almost parallel 
                             columns, and the model associated with the basis 
                             matrix containing those almost parallel rows or 
                             columns.
'''
    
    if (model.IsMIP or model.IsQP or model.IsQCP):
        print("Ill Conditioning explainer only operates on LPs.")
        return None
    if _debug != OFF:
        if _debugger != OFF:
            import pdb; pdb.set_trace()

    if expltype == "COLS":
        expltype = BYCOLS
    elif expltype == "ROWS":
        expltype = BYROWS
    else:
        print("Invalid specification for expltype:", str(expltype) + \
              ": exiting.")
        return None
        
    if relobjtype == "QP":
        relobjtype = SOLVEQP
    elif relobjtype == "LP":
        relobjtype = SOLVELP
    elif relobjtype == "MIP":
        print("MIP subproblem currently unsupported; exiting. ")
        print("Consider setting method to LASSO or submatrix parameter", \
              "to False to obtain smaller explanation.")
        return None
    else:
        print("Invalid specification for relobjtype:", str(relobjtype) + \
              "; exiting.")
        return None

    if method == "DEFAULT":
        method = DEFAULT
    elif method =="ANGLES":
        method = ANGLES
    elif method == "LASSO":
        method = LASSO
    elif method == "RLS":
        method = RLS
    else:
        print("Invalid specification for method:", str(method) + "; exiting.")
        return None
        
        
            

    modvars = model.getVars()      
    modcons = model.getConstrs()   

    if method == ANGLES:
        return angle_explain(model, 1)
#
#   Threshold for ill conditioning is ratio of current feasibility
#   tolerance divided by the current machine precision.  Calculate
#   via logs to avoid dividing small number into larger one.
#
    macheps    = np.finfo(float).eps
    feastol    = model.getParamInfo("FeasibilityTol")[2]
    texp       = math.log(feastol,2) - math.log(macheps, 2)
    condthresh = math.pow(2, texp)     

    splitfreevars = method == LASSO 
    explmodel, splitvardict, RSinginfo, CSinginfo = \
        extract_basis(model, modvars, modcons, expltype, method, \
                      condthresh)
    if explmodel == None:        # Didn't find a basis for original model
        return None              # Nothing to explain
    resmodel  = None
    kappa_stats(model, data, KappaExact)
    if _debug != OFF:
        explmodel.write("explmodel.lp")          
#
#   Minimize the violations of B'y = 0 constraints.  Do not relax
#   the e'y == 1 constraint, and the y variables are free, so they
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
#        v.VBasis = GRB.SUPERBASIC
    explmodel.setParam("NumericFocus", 3)
    explmodel.setParam("Method", 0)
    explmodel.setParam("FeasibilityTol", 1e-9)
    explmodel.setParam("OptimalityTol", 1e-9)
    rowrat = 1.0
    for con in excons[0:nbas]:
        rowvals = []
        rowlhs  = explmodel.getRow(con)
        if rowlhs.size() == 0:
            continue
        rowvals = [rowlhs.getCoeff(j) for j in range(rowlhs.size())]
        freqs = common.get_vector_frequencies(rowvals, 10)
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
            colvals = [col.getCoeff(i) for i in range(col.size())]
            freqs = common.get_vector_frequencies(colvals, 10)  
            rat = freqs[len(freqs) - 1][0] - freqs[0][0] # base 10 exponents
            if rat > colrat:
                colrat = rat
        if colrat > 5:
            explmodel.setParam("ScaleFlag", 2)
    if prmfile != None:
        explmodel.read(prmfile)
    minrelax = method == LASSO or method == RLS
    explmodel.feasRelax(relobjtype, minrelax, None, None, None, \
                        excons[0:nbas], [1.0]*nbas)
    #
    #   Solve configuration completed.  Solve the model that will give us
    #   a certificate of ill conditioning.
    #
    if _debug != OFF:
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
            struct.VBasis = GRB.SUPERBASIC
        for slack in slackvars:
            slack.VBasis = GRB.NONBASIC_LOWER
        cons = explmodel.getConstrs()
        for c in cons:
            c.CBasis = GRB.BASIC
        explmodel.optimize()
        
    #
    #   The y variables were created before calling feasRelax.  Extract
    #   them, as they are the certificate of ill conditioning.
    #
    yvaldict  = {}
    ynamedict = {}
    yvars     = None
    if expltype == BYROWS:
        yvars     = explmodel.getVars()[0:nbas]
        resmodel  = model.copy()
        #
        # Need original model basis statuses as well in order to compute
        # rhs values b - A_Nx_N for the constraints in the explanation
        # model.
        #
        resvars        = resmodel.getVars()
        varsolinfodict = {}
        for ov in model.getVars():
            varsolinfodict[ov.VarName] = (ov.VBasis, ov.X)
        rcons     = resmodel.getConstrs()
        rconsdict = {}
        delcons   = []
        #
        # The order in which we created variables when extracting 
        # the explainer problem typically will not match the order of the
        # constraints in the original model.   So we need to use a dictionary
        # to map the support of the y vector to the correct constraints
        # in the computed explanation.
        #
        combinedcnt = 0
        combinedrhs = 0.0
        combinedlhs = gp.LinExpr()            # y'A; y is inf. certificate
        for c in rcons:
            rconsdict[c.index] = c

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
                suffix_len = len("_GRBPlus")
                yval       = yv.X - splitvardict[yv.index].X
                yname      = yv.VarName[0:-suffix_len]   # original constr name.
            else:
                yval     = yv.X
                yname    = yv.VarName

            zerotol = smalltol
            thiscon  = rconsdict[yv.index]
            rownorm  = L1_rownorm(resmodel, thiscon)
            if rownorm > 0.0:      # Empty lhs is possible
                if smalltol == DEFSMALLTOL:
                    zerotol = min(smalltol, 1/(10*rownorm))
                    if zerotol < macheps:
                        zerotol = 0.0

                
            if abs(yval) <= zerotol:
                delcons.append(rconsdict[yv.index])  # To be filtered out.
            else:
                if _debug == VERBOSE:
                    print("Include constraint ", yname)
                    
                thiscon            = rconsdict[yv.index]
                explname           = "(mult=" + str(yval) + ")" + \
                                     thiscon.ConstrName
                yvaldict[explname] = abs(yval)
                thiscon.ConstrName = explname
                thislhs            = resmodel.getRow(thiscon)
                thisrhs            = thiscon.rhs
                for j in range(thislhs.size()):
                    var     = thislhs.getVar(j)
                    coeff   = thislhs.getCoeff(j)
                    varstat = varsolinfodict[var.VarName][0]
                    if varstat == GRB.BASIC:
                        continue
                    if varstat == GRB.NONBASIC_LOWER:
                        thisrhs  -= (var.LB*coeff)
                    elif varstat == GRB.NONBASIC_UPPER:
                        thisrhs -= (var.UB*coeff)
                    elif varstat == GRB.SUPERBASIC:
                        thisrhs -= (varsolinfodict[var.VarName][1]*coeff)
                    
                thiscon.rhs   = thisrhs
                mult          = yval
                combinedsense = GRB.EQUAL
                # Can be deleted if no issues arise.
                # if thiscon.sense == GRB.GREATER_EQUAL:
                #    mult         = -yval
                #    combinedcnt += 1
                #else:
                #    mult = yval
                #    if thiscon.sense == GRB.LESS_EQUAL:
                #        combinedcnt += 1
                # Inf. Certificate contribution to this constraint and rhs
                combinedlhs.add(thislhs, mult)
                combinedrhs += (thisrhs*mult)

        # if combinedcnt > 0:         can be deleted if no issues arise.
        #    combinedsense = GRB.LESS_EQUAL
        #else:
        #    combinedsense = GRB.EQUAL
        resmodel.addLConstr(combinedlhs, combinedsense, combinedrhs, \
                            COMBINEDROW)  
        resmodel.remove(delcons)
        resmodel.update()
        #
        # Filter out the nonbasic variables of the original model that appear
        # in the constraints ih the explainer model.
        #
        delvars    = []
        resvardict = {}
        resvars    = resmodel.getVars()
        for v in resvars:
            resvardict[v.VarName] = v
        for v in model.getVars():
            if v.vBasis != GRB.BASIC:
                delvars.append(resvardict[v.VarName])
        resmodel.remove(delvars)
        resmodel.update()
        #
        # Print the combined value y'B to the screen, including all nonzero
        # values
        #
        combinedcon = (resmodel.getConstrs())[resmodel.numConstrs - 1]
        print("Vector matrix product of certificate of ill conditioning" + \
              " and basis:")
        print(resmodel.getRow(combinedcon))
        refine_output(resmodel, yvaldict, expltype, submatrix)
    else:              # Column based explanation
        yvars       = explmodel.getVars()[0:nbas]
        resmodel    = explmodel.copy()
        resvars     = resmodel.getVars()
        resvardict   = {}
        delvars     = []
        combinedcol = gp.Column()
        for v in resvars:
            resvardict[v.index] = v
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
                # the original model variable name from the split
                # variable names that have the plus/minus GRB suffixes.
                #
                minusvar   = splitvardict[yv.index]
                yval       = yv.X - minusvar.X
                ynamelen   = len(yv.VarName) - len("_GRBPlus") 
                yname      = yv.VarName[0:ynamelen]   # original variable name
                delvars.append(resvardict[minusvar.index])
            else:
                yval     = yv.X
                yname    = yv.VarName

            thisvar = resvardict[yv.index]
            colnorm = L1_colnorm(resmodel, thisvar)
#
#           By default, with really large row values, we try to capture
#           the significance of even really small row multipliers.
#           But if caller specifies a non default tolerance, we don't
#           overrule them regardless of the presence of really large
#           coefficients that might be missed due to their choice of
#           small multiplier tolerance.
#           
            if smalltol == DEFSMALLTOL:
                zerotol = min(smalltol, 1/colnorm)
                if zerotol < macheps:
                    zerotol = 0.0
            else:
                zerotol = smalltol
            if abs(yval) < zerotol:
                delvars.append(resvardict[yv.index])  # To be filtered out.
            else:
                #
                # Don't include relaxation variables.
                #
                if _debug == VERBOSE:
                    print("Include variable ", yname)
                    
                explname                  = ("(mult=" + str(yval) + ")") + yname
                thisvar.VarName           = explname
                if method == LASSO:  # Bound changed to 0 for Lasso split vars
                    thisvar.LB = -math.inf  
                yvaldict[explname]        = abs(yval)
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
        resmodel.addVar(name=COMBINEDCOL, column=combinedcol)
        rcons = resmodel.getConstrs()
        #
        # e'x = 1 is not in explanation
        # if using lasso approach, then feasrelax phase I constraint
        # also needs to be removed from the explanation
        #
        resmodel.remove(rcons[len(rcons) - 1])
        if splitfreevars:
            resmodel.remove(rcons[len(rcons) - 2])
        resmodel.update()
        refine_output(resmodel, yvaldict, expltype, submatrix)
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
    print("Maximum absolute multiplier value: ", str(max(yvaldict.values()))) 
    print("Minimum absolute multiplier value: ", str(min(yvaldict.values()))) 
    print("--------------------------------------------------------")
    if method == ANGLES:
        return ([], [], None)       # Compatibility with method=ANGLES
    else:
        return resmodel
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
                  method=DEFAULT, condthresh=1e+10):
    #
    #   Does the model have a factorized basis?  If not, need to solve it
    #   first.
    #
    need_basis = False
    try:
        test = model.Kappa
    except AttributeError as e:
        need_basis = True     
    
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
            # TODO: find a better way; this is ugly since it relies on
            # the ability to write to disk.
            #
            model.setParam("IterationLimit", 0)
            model.setParam("Method", 0)
            os.remove("./GRBjunk.bas")
        except gp.GurobiError as e:
            tmp = 0      # No statuses or factorization; solve from scratch
            
        model.optimize() 
        try:       # Confirm basis available after solve.
            tmp2 = modvars[0].VBasis
            tmp1 = modvars[0].X
        except AttributeError as e:    # no basis; no explanation to return
            return None, None, None, None
        
    m        = model.numConstrs
    n        = model.numVars
    #
    #   Extract the basic structural variables from the model into the model
    #   containing the basis matrix of interest.  When explaining by rows,
    #   each column of the basis matrix corresponds to a constraint in the
    #   explainer model, and each row of the basis matrix corresponds to a
    #   variable in the explainer model.
    #   Refinement:  Consider the row and column singletons:
    #   B  = D1   0   0            block row 1 (row singletons)
    #        A1  D2  A2            block row 2 (column singletons)
    #        A3   0  A4            block row 3 (everything else)
    #   Look for special cases where D1 or D2 are individually
    #   ill conditioned (i.e. their ratio of largest to smallest
    #   coefficients exceeds the ill conditioning threshold).
    #   In those cases, just return D1 or D2 as the basis model,
    #   resulting in a much easier computation.
    #
    explmodel   = gp.Model("basismodel")
    RSinginfo    = []            # list of tuples for row singletons
    CSinginfo    = []            # list of tuples for column singletons


    build_explmodel(model, explmodel, modvars, modcons, RSinginfo, CSinginfo, \
                    modeltype)
#
#   Check for special cases where ill conditioning resides just in the
#   the row or column singleton matrices D1 or D2 (i.e. max to min diagonal
#   element exceeds the condition number threshold (~1e+10 for 64 bit doubles)
#
    # TODO: consolidate as much code as possible in this if/elif block.
    if modeltype == BYROWS and len(RSinginfo) > 1:
        maxrat   = 0
        maxcoeff = 0
        mincoeff = float('inf')
        maxvar   = None
        minvar   = None
        maxcon   = None
        mincon   = None
        for diag in RSinginfo:       
            t = abs(diag[COEFF])
            if t > maxcoeff:
                maxcoeff = t
                maxvar = diag[VAR]
                maxcon = diag[CON]
            if t < mincoeff:
                mincoeff = t
                minvar = diag[VAR]
                mincon = diag[CON]
        assert t > 0                   # debug
        if maxcoeff/mincoeff > condthresh:
            #
            # Row singleton matrix D1 is ill conditioned by itself.  Just use
            # the 2x2 matrix associated with the max and min diagonal elements
            # for the explanation
            #
            explmodel.dispose()
            explmodel   = gp.Model("rowsingleton_basismodel")
            build_explmodel(model, explmodel, [minvar, maxvar], \
                           [mincon, maxcon], None, None, BYROWS)
            #
            #
            # The 2x2 matrix is of the form  a1  0 |  0...0
            #                                0  a2 |  0...0
            #                                --------------
            #                                  A1  |   A2   
            # where a1 and a2 are the min or max diagonals.
            # The associated feasrelax model has numerous alterate optimal
            # solutions with minimum relaxation value 1.0.  We seek the one
            # where both rows of the matrix have positive multipliers,
            # so we fix the multipliers to the optimal solution with both
            # multipliers at nonzero values.
            #
            
            denom   = mincoeff + maxcoeff
            minmult = mincoeff/denom
            maxmult = maxcoeff/denom
            explmodel.update()
            for con in explmodel.getConstrs():
                lhs   = explmodel.getRow(con)
                coeff = lhs.getCoeff(0)
                var   = lhs.getVar(0)
                if coeff < maxcoeff:
                    var.LB = maxmult
                    var.UB = maxmult
                else:
                    var.LB = minmult
                    var.UB = minmult
    elif modeltype == BYCOLS and len(CSinginfo) > 1:
        maxrat   = 0
        maxcoeff = 0
        mincoeff = float('inf')
        maxvar   = None
        minvar   = None
        maxcon   = None
        mincon   = None
        for diag in CSinginfo:       
            t = abs(diag[COEFF])
            if t > maxcoeff:
                maxcoeff = t
                maxvar = diag[VAR]
                maxcon = diag[CON]
            if t < mincoeff:
                mincoeff = t
                minvar = diag[VAR]
                mincon = diag[CON]
        assert t > 0                   # debug
        if maxcoeff/mincoeff > condthresh:
            #
            # Column singleton matrix D2 is ill conditioned by itself.  Just use
            # the 2x2 matrix associated with the max and min diagonal elements
            # for the explanation
            #
            explmodel.dispose()
            explmodel   = gp.Model("rowsingleton_basismodel")
            build_explmodel(model, explmodel, [minvar, maxvar], \
                           [mincon, maxcon], None, None, BYROWS)
            #
            # The 2x2 matrix is of the form  a1  0 | A2(1,.)
            #                                0  a2 | A2(2,.) 
            #                                --------------
            #                                0   0  |   A3 
            # where a1 and a2 are the min or max diagonals.
            # The associated feasrelax model has numerous alterate optimal
            # solutions with minimum relaxation value 1.0.  We seek the one
            # where both rows of the matrix have positive multipliers,
            # so we fix the multipliers to the optimal solution with both
            # multipliers at nonzero values.
            #
            denom   = mincoeff + maxcoeff
            minmult = mincoeff/denom
            maxmult = maxcoeff/denom
            explmodel.update()
            for var in explmodel.getVars():
                thiscol = explmodel.getCol(var)
                coeff = thiscol.getCoeff(0)
                if coeff < maxcoeff:
                    var.LB = maxmult
                    var.UB = maxmult
                else:
                    var.LB = minmult
                    var.UB = minmult

#
#   B'y = 0 or By = 0 constraints are done.   All variables are now created, so
#   add the e'y = 1 constraint and we are done.  This is the same,
#   regardless of whether the explanation model is row or column based.
#
    explmodel.update()
    splitvardict = None
    if method == LASSO:
        plusvars     = explmodel.getVars()
        splitvardict = split_mirroredvars(explmodel)
        minusvars    = []
        for v in plusvars:
            minusvar = splitvardict[v.index]
            minusvars.append(minusvar)
#            explmodel.addSOS(GRB.SOS_TYPE1, [v, minusvar]) 
                             
        explmodel.addConstr(gp.quicksum(plusvars) - gp.quicksum(minusvars) == 1)
        explmodel.setObjective(gp.quicksum(plusvars) + gp.quicksum(minusvars))
    else:
        explvars = explmodel.getVars()
        explmodel.addConstr(gp.quicksum(explvars) == 1)
        if method == RLS:
            numvars    = len(explvars)
            sumsquares = gp.QuadExpr()
            sumsquares.addTerms(numvars*[1], explvars, explvars)
            explmodel.setObjective(sumsquares)
                        
                            
            
    explmodel.update()
    return explmodel,splitvardict, RSinginfo, CSinginfo

#
#   Builds the explainer model.  modvars and modcons can be all variables
#   and constraints in the model, but can also be a subset.  Therefore,
#   this routine can be used to build the explainer model after some
#   variables and/or constraints have been filtered out from being
#   possible causes of ill conditioning.
#
def build_explmodel(model, explmodel, modvars, modcons, RSinginfo, CSinginfo,\
                    modeltype):
    starttime        = time.time()
    t1               = 0.0
    t2               = 0.0
    sumtime          = 0.0
    if modeltype == BYROWS:
        #
        #   B'y = 0         
        #   e'y = 1         // normalization of y != 0 constraint
        #   y free
        #
        explvardict = {}
        explcondict = {}
        #
        # Build a dictionary of variables for every constraint in the
        # original model.   Will remove any variables associated with
        # constraints not in modcons after the next loop.
        #
        
        if _debug:
            sumtime  = 0.0
            t1 = time.time()

        conindexdict = {}
        for con in model.getConstrs():
            conindexdict[con.index] = con
            explvardict[con.index]  = explmodel.addVar(lb = -float('inf'), \
                                                       name = con.ConstrName)
        if _debug:
            t2 = time.time()
            sumtime += (t2 - t1)
            print("inner for loop time = ", sumtime)
        
        varlist   = []
        coeflist  = []
        rowcounts = model.NumConstrs * [0]
        for var in modvars:          # structural basic variables
            if var.VBasis != GRB.BASIC:
                continue
            col = model.getCol(var)
            collen   = col.size()
            if collen == 1:          # record column singleton
                if CSinginfo != None:
                    CSinginfo.append((col.getConstr(0), var, col.getCoeff(0)))
            for i in range(collen):
                coeff = col.getCoeff(i)
                con   = col.getConstr(i)
                rowcounts[con.index] += 1
                varlist.append(explvardict[con.index])
                coeflist.append(coeff)
            
            lhs = gp.LinExpr(coeflist, varlist)
#            explcondict[var.VarName] = explmodel.addConstr(lhs == 0, \
#                                                           name=var.VarName)
            explmodel.addConstr(lhs == 0, name=var.VarName)
            coeflist.clear()
            varlist.clear()

        if RSinginfo != None:
            for ind in range(model.NumConstrs):
                if rowcounts[ind] == 1:
                    con = conindexdict[ind]
                    lhs = model.getRow(con)
                    var = lhs.getVar(0)
                    if var.VBasis == GRB.BASIC:
                        RSinginfo.append((con, var, lhs.getCoeff(0)))
            
            
        if CSinginfo != None:
            #
            # slack basic variables.  Ignore if processing row or column
            # singleton model (in which case CSinginfo is None).
            #
            for con in modcons:          # slack basic variables
                if con.CBasis != GRB.BASIC: # TODO:replace next 6 lines with fn.
                    continue
                if con.Sense == GRB.GREATER_EQUAL:
                    coeff = -1.0
                else:
                    coeff = 1.0
                explmodel.addConstr(coeff*explvardict[con.index] == 0, \
                                    name="GRB_slack_" + con.ConstrName)
                #
                # Slacks are column singletons by definition; treat them as
                # part of column singleton diagonal matrix even if also a row
                # singleton.
                #
                CSinginfo.append((con, None, coeff))
        #
        # The previous two loops iterated over modvars from the original
        # model.  If modcons is a subset of the original model constraints,
        # the previous loop may have added variables to the explainer
        # model corresponding to constraints in the original model
        # that are not in the modcons array.  Therefore, remove any
        # variables in the explainer model corresponding to constraints
        # not in the modcons list passed into the routine.
        # passed into the routine.  The explainer model has constraint
        # names based on original model variable names, and variable names
        # based on original model constraint names.
        #
        explmodel.update()
        if len(modcons) < model.NumConstrs:
            modconlist = []
            for con in modcons:
                modconlist.append(con.index)
            delvars = []
            for explvar in explvardict.values():
                if explvar.index in modconlist:
                    continue
                delvars.append(explvar)
            explmodel.remove(delvars)
            explmodel.update()
    else:          # BYCOLS; assumes modeltype has been checked by caller
        #
        #   By = 0        
        #   e'y = 1         // normalization of y != 0 constraint
        #   y free
        #
        explcondict  = {}
        modvarlist   = []
        conindexdict = {}
        for con in model.getConstrs():
            #
            # Constraint initialization for column based explainer model.
            # Build dictionary of constraints for every constraint
            # in the original mode. Will remove any constraints associated
            # with constraints not in modcons after the next loop.
            #
            conindexdict[con.index] = con
            explcondict[con.index] = \
                explmodel.addConstr(0, GRB.EQUAL, 0, name=con.ConstrName)

        #
        # Structural basic variables.  The column from the original model
        # needs to reference the constraint in the explainer model, not
        # the original model.
        #
        colcons   = []
        colcoeffs = []
        rowcounts = model.NumConstrs * [0]
        for var in modvars:          
            if var.VBasis != GRB.BASIC:
                continue
            col     = model.getCol(var)
            collen  = col.size()
            if collen == 1:          # record column singleton
                if CSinginfo != None:
                    CSinginfo.append((col.getConstr(0), var, col.getCoeff(0)))

            for k in range(collen):
                thiscon = col.getConstr(k)
                rowcounts[thiscon.index] += 1
                colcons.append(explcondict[thiscon.index])
                colcoeffs.append(col.getCoeff(k))
            explmodel.addVar(obj=0.0, lb = -float('inf'), name=var.VarName, \
                             column=gp.Column(colcoeffs, colcons))
            colcons.clear()
            colcoeffs.clear()

        if RSinginfo != None:
            for ind in range(model.NumConstrs):
                if rowcounts[ind] == 1:
                    con = conindexdict[ind]
                    lhs = model.getRow(con)
                    var = lhs.getVar(0)
                    if var.VBasis == GRB.BASIC:
                        RSinginfo.append((con, var, lhs.getCoeff(0)))
            
        if CSinginfo != None:
            #
            # slack basic variables.  Ignore if processing row or column
            # singleton model (in which case CSinginfo is None).
            #
            for con in modcons:          # slack basic variables
                if con.CBasis != GRB.BASIC: # TODO: replace next 6 lines with fn.
                    continue
                if con.Sense == GRB.GREATER_EQUAL:
                    coeff = -1.0
                else:
                    coeff = 1.0
                vname = "GRBslack_" + con.ConstrName
                col = gp.Column([coeff], [explcondict[con.index]])
                explmodel.addVar(obj=0.0, lb = -float('inf'), name=vname, \
                                 column=col)
                #
                # Slacks are column singletons by definition; treat them as
                # part of column singleton diagonal matrix even if also a row
                # singleton.
                #
                CSinginfo.append((con, None, coeff))
        #
        # The previous two loops iterated over modvars from the original
        # model.  If modcons is a subset of the original model constraints,
        # the previous loop may have included constraints in the explainer
        # model corresponding to constraint in the original model not in
        # the modcons array.   Therefore, remove any constraints in the
        # explainer model that are not in the modcons list.  The explainer
        # model has variable and constraint names based on the original
        # model variable and constraint names, respectively.
        #
        explmodel.update()
        if len(modcons) < model.NumConstrs:
            modconlist = []
            for con in modcons:
                modconlist.append(con.index)
            delcons = []
            for explcon in explcondict.values():
                if explcon.index in modconlist:
                    continue
                delcons.append(explcon)
            explmodel.remove(delcons)
            explmodel.update()
            
    if _debug != OFF:
        endtime = time.time()
        print("Time build_explmodel = ", endtime - starttime)


def kappa_stats(model, data, KappaExact):
    kappa      = model.Kappa
    kappaexact = -1
    if KappaExact == 1 and (model.numConstrs < 1000000 and \
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
def split_mirroredvars(model, varstosplit=None):
    model.update()
    if varstosplit == None:
        varstosplit = model.getVars()

    newvarlist = []           # This var and dict is used to handle
    newvardict = {}           # quadratic objective and constraints
    for v in varstosplit:
        if v.UB != -v.LB:
            continue
        if v.UB == 0.0 and v.LB == 0.0:       # No need to split vars
            continue                          # fixed at zero
        #
        #  We have a candidate variable to split.
        #
        bnd        = v.UB
        varname    = v.VarName 
        chgvarname = v.VarName  + "_GRBPlus"
        v.VarName  = chgvarname
        v.LB       = 0
        col        = model.getCol(v)
        splitcol   = gp.Column()
        coeflist   = []
        conlist    = []
        for k in range(col.size()):
            coeflist.append(-col.getCoeff(k))
            conlist.append(col.getConstr(k))
        if col.size() > 0:
            splitcol.addTerms(coeflist, conlist)
        newvar = model.addVar(lb=0.0, ub=bnd, obj=-v.obj, vtype=v.Vtype, \
                              name=varname + "_GRBMinus", column=splitcol)
        newvarlist.append(newvar)
        newvardict[v.index] = newvar

    #
    # Linear constraints completed; now update any QCs that contain
    # mirrored variables that need to be split.
    #
    for qcon in model.getQConstrs():
        quadexpr = model.getQCRow(qcon)
        split_quadexpr(quadexpr, newvardict)
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
        split_quadexpr(objexpr, newvardict)
        model.setObjective(objexpr)
    model.update()
    return newvardict
#
#   Takes a quadratic expression and replaces mirrored variables with
#   the appropriate difference of two nonnegative variables.  All
#   mirrored variables are provided in newvardict, which maps the
#   original mirrored variable to its added variable.
    
def split_quadexpr(quadexpr, newvardict):
    for k in range(quadexpr.size()):
        xi = quadexpr.getVar1(k)
        xj = quadexpr.getVar2(k)
        qcoef = quadexpr.getCoeff(k)
        if xi.index in newvardict:      # xi is split into xi(+) - xi(-)
            # xi(-) * xj term
            quadexpr.addTerms(-qcoef, xj, newvardict[xi.index])
            if xj.index in newvardict: # both xi and xj split
                # xi * xj(-) term
                quadexpr.addTerms(-qcoef, xi, newvardict[xj.index])
                # xi(-) * xj(-) term
                quadexpr.addTerms(qcoef, newvardict[xi.index], \
                                      newvardict[xj.index])
        elif xj.index in newvardict: # xj is split, but xi is not.
            # xi * xj(-) term
            quadexpr.addTerms(-qcoef, xi, newvardict[xj.index])


#
#   Look for easy explanations involving almost parallel pairs of
#   matrix rows or columns.   Any pair that is almost parallel can
#   be identified from having a very small angle, i.e., if the
#   difference between the inner product of two vectors and the product
#   if their two norms is within a small tolerance.  Or, if the negative
#   of the inner product and the two norms is within a tolerance, they
#   are almost parallel as well, i.e., the angle is close to 180 degrees.
#   The howmany parameter controls how many near parallel items to report.
#   Default is 1; integer values > 1 specify how many to report; integer
#   values <= 0 indicate report all.  TODO:  This should be callable
#   by the user, so it needs help info.
#
def angle_explain(model, howmany=1, partol=1e-6):
#
#   Help function info  
#    
    '''Specialized ill conditioning explainer for finding pairs of near 
       parallel rows and near parallel columns of the basis matrix.   
       Any basis statuses, in the model will be used, computing the 
       factorization if needed.  If no statuses are available, solve the LP 
       to optimality (or whatever the final solve status is), and use 
       that basis.

       Arguments:
       model      (required) The LP model whose basis will be examined.
                             Basis can be from completed solve or statuses.
                             LP will be optimized if no basis is present.
       howmany    (optional) Number of pairs of near parallel rows and columns
                             to search for.  Default is 1.  Specify a positive
                             integer for a particular number of pairs.   
                             Specify an integer <= 0 to request all pairs.
                             Routine searches rows first, then columns.   
       partol     (optional) Tolerance below which two row or column vectors
                             are considered almost parallel.  Used as a 
                             relative tolerance when comparing the inner
                             product of two vectors with the product of their
                             L1 norms.
       Returns:              returns a list of tuples of the almost parallel 
                             rows and almost parallel columns, and the model 
                             associated with the basis matrix containing those
                             almost parallel rows or columns.
'''

    if _debug != OFF:
        if _debugger != OFF:
            import pdb; pdb.set_trace()
    modcons = model.getConstrs()
    modvars = model.getVars()
    explmodel, splitvardict, junk1, junk2 = extract_basis(model, modvars, \
                                                          modcons, None, False)
    if explmodel == None:             # No basis for original model found.
        return None, None, None       # Nothing to explain
    modcons = explmodel.getConstrs()
    modvars = explmodel.getVars()
    explmodel.remove(modcons.pop(-1))      # Just want basis matrix
    explmodel.update()         
    #
    # Calculate mod value for hashing function, which depends on
    # max possible number of bits for an integer (typically 64).
    # Also set up some dictionaries for hashing calculations.
    #
    modval      = 1 + math.ceil(math.log(sys.maxsize,2))
    rowhashdict = {}
    for con in modcons:
        rowhashdict[con.ConstrName] = 0
    colhashdict = {}
    for var in modvars:
        colhashdict[var.VarName] = 0

    #
    # Look for almost parallel rows.   Hash the rows as follows.
    # For each row define a 64 bit int h[i] such that if variable
    # xj appears in the row with a coefficient larger than epsilon
    # you set bit j % 64 in h_i. Then, for a pair (i1,i2) of rows,
    # if h[i1] != h[i2], the two rows cannot be almost parallel.
    # First, create the hash values and put the constraints with the
    # same hash value in the same bucket.   Don't hash small nonzeros
    # below partol; that way we pick up the case where one row is
    # (1 1) and another intersects the same two nonzeros but is
    # (1 1 eps).   We want to flag those as almost parallel, so we
    # want them to have the same hash value.  We do this in a normalized
    # fashion, so that we will get the same conclusion for (1 1), (1, 1 eps)
    # and (1000 1000), (1000, 1000, 1000*eps).  Similarly, we don't
    # consider rows of the form (1e-5, 1e-5, 0) and (1e-5, 1e-5, 1e-7)
    # almost parallel since they are significantly different in a relative
    # sense.
    #
    bucketlists = [[] for i in range(modval)]
    for con in modcons:
        rowlhs  = explmodel.getRow(con)
        normvec = []
        for j in range(rowlhs.size()):
            normvec.append(rowlhs.getCoeff(j))
        rownorm = L1_norm(normvec)
                       
        for j in range(rowlhs.size()):
            var  = rowlhs.getVar(j)
            coef = rowlhs.getCoeff(j)
            if abs(coef) > partol*rownorm:
                varbit = var.index % (modval - 1)
                rowhashdict[con.ConstrName] |= (1 << varbit)
        bucketlists[rowhashdict[con.ConstrName] % (modval - 1)].append(con)
    #
    # Check the matrix rows in the same bucket to see if they are
    # almost parallel.
    #
    almostparcount   = 0
    almostparrowlist = []
    quittingtime     = False
    starttime        = time.time()    
    for bucket in bucketlists:
        for k in range(len(bucket)):
            lhs1      = explmodel.getRow(bucket[k])
            size1     = lhs1.size()     
            norm1     = 0.0
            row1array = explmodel.numVars*[0.0]
            #
            # Treat first row as a dense vector, second row as a sparse one.
            #
            for i in range(size1):
                ind1            = lhs1.getVar(i).index
                coef1           = lhs1.getCoeff(i) 
                row1array[ind1] = coef1
                norm1          += coef1**2
            for j in range((k+1),len(bucket)):
                lhs2  = explmodel.getRow(bucket[j])   
                size2 = lhs2.size()
                if size1 != size2:
                    continue               
                #
                # Simple length based filtering failed.  Compute inner products
                # and norms.  Need to handle the case where two rows intersect
                # the same variables, but in a different order.
                #
                #
                # Compute the angle.  TODO: Make this calculation more
                # precise.  Consider Kahan summation, and avoid adding
                # big numbers and much smaller ones.
                #
                dotprod  = 0.0
                norm2    = 0.0
                skip     = False
                for i in range(size2):
                    ind2 = lhs2.getVar(i).index
                    if row1array[ind2] != 0.0:
                        coef2    = lhs2.getCoeff(i)
                        dotprod += row1array[ind2]*coef2
                        norm2   += coef2**2
                    else:  # Intersecting variables differ; nothing to see here
                        skip=True
                        break
                
                if skip:  # Early exit due to index or coeff mismatch
                    continue
                norm1 = math.sqrt(norm1)    
                norm2 = math.sqrt(norm2)
                #
                # Normalize the tolerance in the test; want the same
                # result for vectors u and v and 1000u and 1000v
                #
                if abs(abs(dotprod) - norm1*norm2) > partol*norm1:
                    continue
                #
                # Found two almost parallel rows.
                #
                almostparcount += 1
                almostparrowlist.append((bucket[k], bucket[j]))
                if howmany <= 0 or almostparcount < howmany:
                    continue
                quittingtime=True
                break                       # Exit for j loop
            if quittingtime:
                break                       # Exit for k loop
        if quittingtime:
            break                           # Exit for bucket loop

    if _debug != OFF:    
        endtime = time.time()
        print("Time for main row based loop = ", endtime - starttime)
    del row1array
    #
    # Now look for almost parallel columns.   Hash the rows in a smaller
    # manner, except this time the bit values are based on row indices
    # rather than column indices.
    #
    bucketlists = [[] for i in range(modval)]
    for var in modvars:
        col = explmodel.getCol(var)
        normvec = []
        for i in range(col.size()):
            normvec.append(col.getCoeff(i))
        colnorm = L1_norm(normvec)
        for i in range(col.size()):
            con  = col.getConstr(i)
            coef = col.getCoeff(i)
            if abs(coef) > partol*colnorm:
                varbit = con.index % (modval - 1)
                colhashdict[var.VarName] |= (1 << varbit)
        bucketlists[colhashdict[var.VarName] % (modval - 1)].append(var)
    #
    # Check the matrix columns in the same bucket to see if they are
    # almost parallel.
    #
    almostparcollist = []
    quittingtime     = False
    starttime        = time.time()
    for bucket in bucketlists:
        for k in range(len(bucket)):
            col1  = explmodel.getCol(bucket[k])
            col1array = explmodel.numVars*[0.0]
            norm1     = 0.0
            size1 = col1.size()     
            for i in range(size1):    # Treat first column as dense vector.
                ind1            = col1.getConstr(i).index
                coef1           = col1.getCoeff(i)
                col1array[ind1] = coef1
                norm1          += coef1**2
            for j in range((k+1),len(bucket)):
                col2  = explmodel.getCol(bucket[j])   
                size2 = col2.size()
                if size1 != size2:
                    continue               
                #
                # Simple length based filtering failed.  Compute inner products
                # and norms.  Need to handle the case where two columns
                # intersect the same rows, but in a different order.
                #
                # Compute the angle.  TODO: Make this calculation more
                # precise.  Consider Kahan summation, and avoid adding
                # big numbers and much smaller ones.
                #
                dotprod  = 0.0
                norm2    = 0.0
                skip     = False
                for i in range(size2):
                    ind2  = col2.getConstr(i).index
                    if col1array[ind2] != 0.0:
                        coef2    = col2.getCoeff(i)
                        dotprod += col1array[ind2] * coef2
                        norm2   += coef2**2
                    else:
                        #
                        # Intersecting constraints differ; move along, nothing
                        # to see here.
                        #
                        skip=True
                        break
                
                if skip:  # Early exit due to index or coeff mismatch
                    continue
                norm1 = math.sqrt(norm1)    
                norm2 = math.sqrt(norm2)
                #
                # Normalize the tolerance in the test; want the same
                # result for vectors u and v and 1000u and 1000v
                #
                if abs(abs(dotprod) - norm1*norm2) > partol*norm1:
                    continue
                #
                # Found two almost parallel rows.
                #
                almostparcount += 1
                almostparcollist.append((bucket[k], bucket[j]))
                if howmany <= 0 or almostparcount < howmany:
                    continue
                quittingtime=True
                break                       # Exit for j loop
            if quittingtime:
                break                       # Exit for k loop
        if quittingtime:
            break                           # Exit for bucket loop

    if _debug != OFF:
        endtime = time.time()
        print("Time for main column based loop = ", endtime - starttime)
    
    return almostparrowlist, almostparcollist, explmodel




#                                     
#   Improve the readability of the output
#   1) List the constraints in the explanation sorted by largest
#      absolute row multiplier first.
#   2) If model has row or column singletons, include only selected
#      submatrices in the explanations (details in comments below).
#   
def refine_output(model, absmultdict, expltype, submatrix):
#
#   model is the explainer model, not the original one.
#   This routine modifies the model, with the constraints sorted in
#   descending order of the values in the absmultdict dictionary.
#
#    import pdb; pdb.set_trace()
    quit = False
    if expltype == BYROWS:
        condict       = {}
        #
        # Skip the combined constraint (the last one); it has no multiplier
        #
        cons          = model.getConstrs()[0:model.numConstrs - 1]
        concnt        = 0
        for c in cons:
            condict[c.ConstrName] = c
            concnt += 1
        #
        # Sort the dictionary by absolute row multiplier value from
        # highest to lowest.  Dictionaries cannot be sorted directly,
        # so we use the sorted function to obtain a list of tuples
        # sorted by the dictionary value, then transform that list of
        # tuples back into a dictionary.  Details are at
        # https://www.freecodecamp.org/news/sort-dictionary-by-value-in-python/
        #
        sorteddict = dict(sorted(absmultdict.items(), key=lambda x:x[1], \
                                 reverse=True))
        #
        # Add the constraints in sorted order to the end of the model, then
        # remove their duplicates in unsorted order that are the first concnt
        # constraints
        #
        for conname in sorteddict.keys():
            contoadd = condict[conname]
            lhs      = model.getRow(contoadd)
            model.addConstr(lhs, contoadd.sense, contoadd.rhs, \
                            name=contoadd.ConstrName)
        model.remove(cons[0:concnt])
        model.update()
    elif expltype == BYCOLS:
        vardict       = {}
        #
        # Skip the combined column (the last one); it has no multiplier
        #
        vars          = model.getVars()[0:model.numVars - 1]
        varcnt        = len(vars)
        for v in vars:
            vardict[v.VarName] = v
        #
        # Sort the dictionary by absolute col multiplier value from
        # highest to lowest.   Analogous to row multiplier sorting described
        # in detail in the comment in the BYROWS code above
        #
        sorteddict = dict(sorted(absmultdict.items(), key=lambda x:x[1], \
                                 reverse=True))
        #
        # Add the variables in sorted order to the end of the model, then
        # remove their duplicates in unsorted order that are the first varcnt
        # variables
        #
        for varname in sorteddict.keys():
            vartoadd = vardict[varname]
            model.addVar(lb=vartoadd.lb, ub=vartoadd.ub, \
                         column=model.getCol(vartoadd), name=vartoadd.VarName)
        model.remove(vars[0:varcnt])
        model.update()
    #
    # First refinement complete.  Model constraints or variables will
    # now be listed by largest absolute multiplier value first.
    #
    
    if submatrix:
             
        #
        # Additional reductions in explanation that may permit it to
        # be an ill conditioned submatrix of the original basis matrix.
        # Details depend on whether explanation is row or column based.
        #
        #   B  = D1   0   0            block row 1 (row singletons)
        #        A1  D2  A2            block row 2 (column singletons)
        #        A3   0  A4            block row 3 (everything else)
        #
        # Note that if we arrive here, preprocessing has already determined
        # that D1 and D2 individually are not ill conditioned. 
        #
        # First, perform tasks common to both row and column based refinements.
        # Build lists of variables and constraints in the various
        # blocks.  Original singleton data structures RSinginfo and CSinginfo
        # were computed during basis extraction, so we need to remove the
        # associated constraints and variables that are not part of the
        # actual explanation.  Easier to just rebuild these data structures
        # than use the original ones.
        #
        cons       = model.getConstrs()
        vars       = model.getVars()
        D1condict  = {}
        D1vardict  = {}
        D2condict  = {}
        D2vardict  = {}
        D1coefdict = {}
        D2coefdict = {}
        for c in cons:
            lhs = model.getRow(c)
            if lhs.size() != 1:
                continue
            #
            # Row singleton found
            #
            var = lhs.getVar(0)
            if var.VarName.startswith("GRB_slack"):
                #
                # Treat row singleton slacks as column singletons
                # rather than row singletons.
                #
                D2condict[c.ConstrName] = c
                D2vardict[var.VarName]  = None
                D2coefdict[var.VarName] = 1.0      # sign doesn't matter
            else:               # Row singleton, part of D1
                D1condict[c.ConstrName]  = c
                D1vardict[var.VarName]   = var
                D1coefdict[c.ConstrName] = lhs.getCoeff(0)
                
        for var in vars:
            if var.VarName.startswith("GRB_slack"):
                continue           # Slack column singletons already processed.
            col    = model.getCol(var)
            if col.size() != 1:
                continue
            #
            # Record column singleton info
            #
            con = col.getConstr(0)
            D2condict[con.ConstrName] = con
            D2condict[var.VarName]    = var
            D2coefdict[var.VarName]   = col.getCoeff(0)

        if expltype == BYROWS:
            #
            #   B  = D1   0   0            block row 1 (row singletons)
            #        A1  D2  A2            block row 2 (column singletons)
            #        A3   0  A4            block row 3 (everything else)
            #
            # Reduction #1:  Let w1, w2 and w3 be the row multipliers from the
            # certificate of ill conditioning associated with blocks 1, 2 and
            # 3 respectively.   If w3'A4 shows that A4 is ill conditioned, then
            # we know that B is ill conditioned by setting w2 = 0,
            # w1 = -w3'D1^(-1)A3.   So just output the rows and columns in A4 by
            # removing the variables and constraints associated with blocks 1
            # and 2, and use the remaing  submatrix A4 as the explanation.
            #
            combocon = None
            for c in cons:
                if c.ConstrName == COMBINEDROW:
                    combocon = c
                    break
            combolhs      = model.getRow(combocon)
            A4combovars   = []
            A4combocoeffs = []
            for j in range(combolhs.size()):
                var   = combolhs.getVar(j)
                if var.VarName in D1vardict or var.VarName in D2vardict:
                    continue
                #
                # This term of the combined constraint is a column of A4
                #
                A4combovars.append(var)
                A4combocoeffs.append(abs(combolhs.getCoeff(j))) # contains w3'A4
            #
            # We have w3'A4.  We need w3 so we can assess near singularity 
            # of A4. But make sure A4 has positive dimension before doing so.
            #
            if len(D1condict) + len(D2condict) < model.NumConstrs:
                w3vals = []
                for c in cons:
                    if c.ConstrName in D1condict or c.ConstrName in D2condict:
                        continue
                    if c.ConstrName.startswith(COMBINEDROW):
                        continue
                    w3vals.append(absmultdict[c.ConstrName])
                w3norm      = L1_norm(w3vals)
                A4combonorm = L1_norm(A4combocoeffs)
                if A4combonorm/w3norm < 0.01:
                    #
                    # A4 ill conditioned.  Remove the constraints in block rows
                    # 1 and 2, and the variables in D1 and D2 from the explainer
                    # model.
                    #
                    delcons = list(D1condict.values())
                    delcons.append(list(D2condict.values()))
                    delvars = list(D1vardict.values())
                    delvars.append(list(D2vardict.values()))
                    model.remove(delcons)
                    model.update()
                    model.remove(delvars)
                    model.update()
                    #
                    # The remaining explanation model is just the variables
                    # and constraints in A4.   Don't look for additional
                    # refinements.
                    #
                    quit = True
        elif expltype == BYCOLS:
            #
            #   B  = D1   0   0            block row 1 (row singletons)
            #        A1  D2  A2            block row 2 (column singletons)
            #        A3   0  A4            block row 3 (everything else)
            #
            #
            # Reduction #1:  Let w1, w2 and w3 be the column multipliers
            # from the certificate of ill conditioning associated with 
            # blocks 1, 2 and 3 respectively.   If A4w3 shows that A4
            # is ill conditioned, then we know that B is ill conditioned
            # by setting w1 = 0, w2 = -D2^(-1)A2. So just output the rows
            # and columns in A4 by removing the variables and constraints
            # associated with blocks 1 and 2, and use the remaining
            # submatrix A4 as the explanation.
            #
            combovar = None
            for v in vars:
                if v.VarName == COMBINEDCOL:
                    combovar = v
                    break;
            combocol      = model.getCol(combovar)
            A4combocons   = []
            A4combocoeffs = []
            for i in range(combocol.size()):
                con = combocol.getConstr(i)
                if con.ConstrName in D1condict or con.ConstrName in D2condict:
                    continue
                #
                # This term of the dombined column is a row of A4
                #
                A4combocons.append(con)
                A4combocoeffs.append(abs(combocol.getCoeff(i)))
            #
            # We have A4w3.  We need w3 so we can assess near singularity of
            # A4.   But make sure  A4 has positive dimension before doing so.
            #
            if len(D1vardict) + len(D2vardict) < model.NumVars:
                w3vals = []
                for v in vars:
                    if v.VarName in D1vardict or v.VarName in D2vardict:
                        continue
                    if v.VarName.startswith(COMBINEDCOL):
                        continue
                    w3vals.append(absmultdict[v.VarName])
                w3norm      = L1_norm(w3vals)
                A4combonorm = L1_norm(A4combocoeffs)
                if A4combonorm/w3norm < 0.01:
                    #
                    # A4 ill conditioned.  Remove the constraints and variables
                    # associated with blocks 1 and 2.
                    # TODO: If this stays indentical to the corresponding
                    # code in the BYROWS refinement, put this block into
                    # a single function.
                    #
                    delcons = list(D1condict.values())
                    delcons.append(list(D2condict.values()))
                    delvars = list(D1vardict.values())
                    delvars.append(list(D2vardict.values()))
                    model.remove(delcons)
                    model.update()
                    model.remove(delvars)
                    model.update()
                    #
                    # The remaining explanation model is just the variables
                    # and constraints in A4.   Don't look for additional
                    # refinements.
                    #
                    quit = True

                                  
                
            
    if not quit:
        #
        # Preprocessing checked whether D1 or D2 individually was ill
        # conditioned and took care of the explanation.   Here we check
        # if D1 and D2 combined are ill conditioned even though individually
        # they are well conditioned (e.g., D1 has all values around 1e-5
        # while D2 has all values around 1e+5).
        #
        quit = True
        
            
    
    
            
#                                     
#   Improve the readability of the output.
#   1) List the variables in the explanation sorted by largest
#      absolute column multiplier first. 
#   
def refine_col_output(model, abscolmultdict):
#
#   model is the explainer model, not the original one.
#   This routine modifies the model, with the variables sorted in
#   descending order of the values in the abscolmultdict dictionary.
#
#    import pdb; pdb.set_trace()
    vardict       = {}
    #
    # Skip the combined column (the last one); it has no multiplier
    #
    vars          = model.getVars()[0:model.numVars - 1]
    varcnt        = len(vars)
    for v in vars:
        vardict[v.VarName] = v
    #
    # Sort the dictionary by absolute col multiplier value from
    # highest to lowest.  Dictionaries cannot be sorted directly,
    # so we use the sorted function to obtain a list of tuples
    # sorted by the dictionary value, then transform that list of
    # tuples back into a dictionary.  Details are at
    # https://www.freecodecamp.org/news/sort-dictionary-by-value-in-python/
    #
    sorteddict = dict(sorted(abscolmultdict.items(), key=lambda x:x[1], \
                             reverse=True))
    #
    # Add the variables in sorted order to the end of the model, then
    # remove their duplicates in unsorted order that are the first varcnt
    # variables
    #
    for varname in sorteddict.keys():
        vartoadd = vardict[varname]
        model.addVar(lb=vartoadd.lb, ub=vartoadd.ub, \
                     column=model.getCol(vartoadd), name=vartoadd.VarName)
    model.remove(vars[0:varcnt])
    model.update()
    #
    # Model variables will now be listed by largest absolute multiplier
    # value first.
    #
    

#
#   L1 norm of a specific constraint in a model.
#    
def L1_rownorm(model, con):
    lhs = model.getRow(con)
    normvec = []
    for j in range(lhs.size()):
        normvec.append(lhs.getCoeff(j))
    result = L1_norm(normvec)
    return result

#
#   L1 norm of a specific variable in a model.
#    
def L1_colnorm(model, var):
    col = model.getCol(var)
    normvec = []
    for j in range(col.size()):
        normvec.append(col.getCoeff(j))
    result = L1_norm(normvec)
    return result    
#
#   L1 norm of an arbitrary vector.  Don't need super accurate
#   values here, so don't worry about Kahan summation or other
#   tactics for precision
#
def L1_norm(vec):
    sum = 0.0
    for v in vec:
        sum += abs(v)
    return sum
#
#   TODO: include Skeel condition number calculation
#
    

#
#   Debug routine
#
def bucket_info(bucketlists):
    histogram = []
    for bucket in bucketlists:
        histogram.append(len(bucket))
    return histogram
