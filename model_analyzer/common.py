import gurobipy as gp
from collections import defaultdict
import math


def get_vtype(var):
    if var.VType == gp.GRB.INTEGER and var.LB >= 0 and var.UB <= 1:
        return gp.GRB.BINARY
    else:
        return var.VType


def max_value(variable, coefficient):
    """
    Calculates the maximum attainable value for a variable. Either the upper bound or infinity.
    :param variable: The variable
    :param coefficient: The coefficient with which the upper bound should multiplied
    :return: The maximum value
    """
    if coefficient == 0:
        return 0
    elif coefficient > 0:
        return min(gp.GRB.INFINITY, variable.UB * coefficient)
    else:
        return min(gp.GRB.INFINITY, variable.LB * coefficient)


def min_value(variable, coefficient):
    """
    Calculates the minimum attainable value for a variable. Either the upper bound or infinity.
    :param variable: The variable
    :param coefficient: The coefficient with which the lower bound should multiplied
    :return: The minimum value
    """
    if coefficient == 0:
        return 0
    elif coefficient > 0:
        return max(-gp.GRB.INFINITY, variable.LB * coefficient)
    else:
        return max(-gp.GRB.INFINITY, variable.UB * coefficient)


def numfinite(vals):
    """
    Counts the number of elements in an array take on finite values
    :param vals:  The array of values to check
    :return: True or False
    """
    count = 0
    for v in vals:
        if v > -gp.GRB.INFINITY and v < gp.GRB.INFINITY:
            count += 1
    return count


def isfinite(vals):
    """
    Checks if all elements in an array take on finite values
    by calling the numfinite routine to count the number of finite
    values.
    :param vals:  The array of values to check
    :return: True or False
    """

    return numfinite(vals) == len(vals)


def qterm_min_value(var1, var2, qcoeff):
    """
    Calculates the minimum attainable value for a single quadratic term
    More complicated than the corresponding min_value routine for linear
    terms, as we need to consider 4 lower/upper bound combinations for the
    two variables in the quadratic term.
    :param var1, var2: The two variables comprising the quadratic term
    :param qcoeff: The coefficient of the quadratic term
    :return: The minimum value of the quadratic term
    """
    if qcoeff == 0:
        return 0
    #   TODO: Handle the case of infeasible bounds where at least one
    #   lower bound exceeds the corresponding upper bound, regardless of
    #   violation amount.
    #    import pdb; pdb.set_trace()

    l1 = var1.LB
    l2 = var2.LB
    u1 = var1.UB
    u2 = var2.UB
    vallist = []

    """
    Start with the easy cases that are similar to linear terms, i.e.
    none of the 4 bounds have opposite signs.
    """

    if l1 >= 0 and l2 >= 0:
        if u1 >= gp.GRB.INFINITY or u2 >= gp.GRB.INFINITY:
            vallist.extend([l1 * l2, gp.GRB.INFINITY])
        else:  # all 4 bounds are finite
            vallist.extend([l1 * u2, l2 * u1, l1 * l2, u1 * u2])
    elif u1 <= 0 and u2 <= 0:
        if l1 <= -gp.GRB.INFINITY or l2 <= -gp.GRB.INFINITY:
            vallist.extend([u1 * u2, gp.GRB.INFINITY])
        else:  # all 4 bounds are finite
            vallist.extend([l1 * u2, l2 * u1, l1 * l2, u1 * u2])

    else:
        """
        Now consider the harder cases that are probably rare in practice,
        where not all 4 variable bounds are nonnegative or nonpositive.
        Need to consider all 4 possible variable bound products, and handle
        cases of one or more infinite bounds separately.
        """
        bndlist = [l1, l2, u1, u2]
        lblist = [l1, l2]
        ublist = [u1, u2]
        minval = gp.GRB.INFINITY
        maxval = -gp.GRB.INFINITY
        if isfinite(bndlist):  # Easy; all bilinear products finite
            vallist = [l1 * u2, l2 * u1, l1 * l2, u1 * u2]
        else:
            """
            At least one bound at +-infinity.  Rather than dealing with
            all the possible different cases, there really are only a handful
            of cases to consider that we can categorize based mostly on the
            number of infinite bounds among the list of 4 bounds for var1
            and var2
            """
            ninfbds = len(bndlist) - numfinite(bndlist)
            if ninfbds >= 3:
                """
                At least one xi*xj term at +inf and at least one at -inf,
                so we know the max and min bound products without much work
                """
                vallist.extend([-gp.GRB.INFINITY, gp.GRB.INFINITY])
            elif ninfbds == 2:
                """
                We have one special case here.   If var1 or var2 are free,
                then at least one xi*xj term is +inf and at least one is
                -inf regardless of the signs of the lower and upper bound
                of the other variable.
                """
                if numfinite([l1, u1]) == 0 or numfinite([l2, u2]) == 0:
                    vallist.extend([-gp.GRB.INFINITY, gp.GRB.INFINITY])
                elif numfinite([l1, l2]) == 0:
                    """
                    Both lower bounds negative infinite, so both upper bounds
                    finite.   Min and Max xi*xj terms depend on signs of
                    the finite bounds.
                    """
                    if u1 >= 0 or u2 >= 0:
                        vallist.extend([-gp.GRB.INFINITY, gp.GRB.INFINITY])
                    else:  # u1, u2 nonpositive; their product is the min xi*xj
                        vallist.extend([u1 * u2, gp.GRB.INFINITY])
                elif numfinite([u1, u2]) == 0:
                    """
                    Both upper bounds infinite, so both lower bounds
                    finite.   Min and Max xi*xj terms depend on signs of
                    the finite bounds.
                    """
                    if l1 <= 0 or l2 <= 0:
                        vallist.extend([-gp.GRB.INFINITY, gp.GRB.INFINITY])
                    else:  # l1, l2 positive; their product is the min xi*xj
                        vallist.extend([l1 * l2, gp.GRB.INFINITY])
                else:
                    """
                    One lower bound is -infinite and one upper bound is
                    infinite for the other variable (we handled the
                    free variable case previously).  If finite bounds are
                    both nonnegative or nonpositive, then min and max
                    xi*xj values are -inf and +inf.   Otherwise, need to
                    consider different locations of finite bounds of
                    opposite sign.
                    """
                    if numfinite([l1, u2]) == 2:
                        if l1 <= 0 and u2 <= 0 or l1 >= 0 and u2 >= 0:
                            vallist.extend([-gp.GRB.INFINITY, gp.GRB.INFINITY])
                        else:
                            if u1 >= 0 and u2 >= 0 and l1 <= 0 and l2 <= 0:
                                vallist.extend([-gp.GRB.INFINITY, gp.GRB.INFINITY])
                            else:
                                vallist.extend([-gp.GRB.INFINITY, l1 * u2])
                    else:  # l2 and u1 finite
                        if l2 <= 0 and u1 <= 0 or l2 >= 0 and u1 >= 0:
                            vallist.extend([-gp.GRB.INFINITY, gp.GRB.INFINITY])
                        else:
                            if u1 >= 0 and u2 >= 0 and l1 <= 0 and l2 <= 0:
                                vallist.extend([-gp.GRB.INFINITY, gp.GRB.INFINITY])
                            else:  # l2 >= 0 and u1 <= 0
                                vallist.extend([-gp.GRB.INFINITY, l2 * u1])
            else:
                """
                All that remains is the case when one bound is +-infinity. We've already
                handled the two cases of l1 and l2 nonnegative and u1 and u2 nonpositive.
                Find the +- infinite bound, then look at the signs of the two bounds that
                have products with the infinite bound to determine the possible products
                of pairs of variable bounds
                """
                if numfinite([l1, u1]) == 2:  # x2 has the +- infinite bound
                    if l1 <= 0 and u1 >= 0:
                        vallist.extend([-gp.GRB.INFINITY, gp.GRB.INFINITY])
                    else:  # l1 and u1 have the same sign
                        infval = math.copysign(1, u1) * gp.GRB.INFINITY
                        if isfinite([l2]):  # u2 infinite
                            vallist.extend([l1 * l2, l2 * u1, infval])
                        else:  # l2 -infinite
                            vallist.extend([l1 * u2, u1 * u2, -infval])

                else:  # x1 has the +- infinite bound
                    if l2 <= 0 and u2 >= 0:
                        vallist.extend([-gp.GRB.INFINITY, gp.GRB.INFINITY])
                    else:  # l2 and u2 have the same sign
                        infval = math.copysign(1, u2) * gp.GRB.INFINITY
                        if isfinite([l1]):  # u1 infinite
                            vallist.extend([l1 * l2, l1 * u2, infval])
                        else:  # l1 -infinite
                            vallist.extend([u1 * u2, l2 * u1, -infval])

    """
    The heavy lifting is done.   Just need to look at min and max
    values of the at most 4 possible different products of the two variable bounds
    to determine min quadratic term value
    """
    minval = min(vallist)
    maxval = max(vallist)
    if qcoeff > 0:
        if isfinite([minval]):
            return qcoeff * minval
        else:
            return minval
    else:
        if isfinite([maxval]):
            return qcoeff * maxval
        else:
            return -maxval


# TODO: Need a qterm_max_value.


def utc_date_formatted(dt):
    """
    Format the date into a string
    :param dt: The date
    :return: The string
    """
    weekday = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"][dt.weekday()]
    month = [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ][dt.month - 1]
    return f"{weekday}, {dt.day:02d} {month} {dt.year:04d} {dt.hour:02d}:{dt.minute:02d}:{dt.second:02d} GMT"


def count_decimal_digits(number):
    """
    Count the number of decimal digits in a number
    :param number: The number
    :return: The number of decimal digits
    """
    number_parts = repr(number).split(".")
    if len(number_parts) == 1:
        return 0
    else:
        if number_parts[1] == "0":
            return 0
        else:
            return len(number_parts[1].rstrip("0"))


def get_variable_bounds(m, data, basis):
    coeff_lb_count = defaultdict(int)
    coeff_ub_count = defaultdict(int)

    min_exponent = float("inf")
    max_exponent = -float("inf")

    for var in m.getVars():
        if var.LB != 0 and abs(var.LB) < gp.GRB.INFINITY:
            exponent = int(math.floor(math.log(abs(var.LB), basis)))
            coeff_lb_count[exponent] += 1
            if exponent < min_exponent:
                min_exponent = exponent
            if exponent > max_exponent:
                max_exponent = exponent

        if var.UB != 0 and abs(var.UB) < gp.GRB.INFINITY:
            exponent = int(math.floor(math.log(abs(var.UB), basis)))
            coeff_ub_count[exponent] += 1
            if exponent < min_exponent:
                min_exponent = exponent
            if exponent > max_exponent:
                max_exponent = exponent

    if min_exponent == float("inf"):
        return [], []

    else:
        coeff_lb_result = [
            [exponent, coeff_lb_count[exponent]]
            for exponent in range(min_exponent, max_exponent + 1)
        ]
        coeff_ub_result = [
            [exponent, coeff_ub_count[exponent]]
            for exponent in range(min_exponent, max_exponent + 1)
        ]
        return coeff_lb_result, coeff_ub_result


def get_vector_frequencies(vec, base=10):
    count = defaultdict(int)
    result = []

    for v in vec:
        if v != 0 and abs(v) != float("inf"):
            count[int(math.floor(math.log(abs(v), base)))] += 1

    if len(count) > 0:
        result = [
            [exponent, count[exponent]]
            for exponent in range(min(count.keys()), max(count.keys()) + 1)
        ]

    return result


def get_obj_c_frequencies(model, base=10):
    int_result = []
    cont_result = []
    intobjcoeffs = []
    contobjcoeffs = []
    for x in model.getVars():
        if x.Obj != 0:
            if x.VType == "C":
                contobjcoeffs.append(x.Obj)
            else:
                intobjcoeffs.append(x.Obj)

    if len(contobjcoeffs) > 0:
        cont_result = get_vector_frequencies(contobjcoeffs, base)
    if len(intobjcoeffs) > 0:
        int_result = get_vector_frequencies(intobjcoeffs, base)

    return cont_result, int_result


#
#   Extract linear constraints with names specified by prefix and optional
#   suffix
#


def extractconstraints(model, prefix, suffix=None):
    cons = []
    for c in model.getConstrs():
        cstr = c.ConstrName
        if prefix != None:
            if not cstr.startswith(prefix):
                continue
        if suffix != None:
            if not cstr.endswith(suffix):
                continue
        cons.append(c)

    return cons


#
#   Extract quadratic constraints with names specified by prefix
#   and optional suffix
#


def extractqconstraints(model, prefix, suffix=None):
    qcons = []
    for qc in model.getQConstrs():
        cstr = qc.QCName
        if prefix != None:
            if not cstr.startswith(prefix):
                continue
        if suffix != None:
            if not cstr.endswith(suffix):
                continue
        qcons.append(qc)

    return qcons


#
#   Extract variables with names specified by prefix and optional suffix
#   Also allow for extracting by upper or lower bounds above/below a
#   specified threshold.  Also allow for extracting by a string of types
#   consisting of a subset of "BIC" (Binary/Integer/Continuous)
#


def extractvariables(model, prefix=None, suffix=None, bndthresh=None, types=None):
    vars = []
    for v in model.getVars():
        vstr = v.VarName
        if prefix != None:
            if not vstr.startswith(prefix):
                continue
        if suffix != None:
            if not vstr.endswith(suffix):
                continue
        if bndthresh != None:
            if (v.UB >= GRB.INFINITY or v.UB < bndthresh) and (
                v.LB <= -GRB.INFINITY or v.LB > -bndthresh
            ):
                continue
        if types != None:
            hit = False
            for char in types:
                if v.Vtype == char:
                    hit = True
                    break
            if not hit:
                continue
        vars.append(v)

    return vars


#
#   Returns the variables contained in a linear expression
#
def extractvariablesfromexpr(linexpr):
    retvars = []
    for j in range(linexpr.size()):
        retvars.append(linexpr.getVar(j))

    return retvars


#
#   Iterate through the dictionary in data, which is assumed
#   to contain string keys and  integer values, and print the results
#
def printvals(datadict):
    for datakey in datadict:
        print(datakey, " = ", str(datadict[datakey]))
