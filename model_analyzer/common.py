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


def utc_date_formatted(dt):
    """
    Format the date into a string
    :param dt: The date
    :return: The string
    """
    weekday = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"][dt.weekday()]
    month = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
             "Oct", "Nov", "Dec"][dt.month - 1]
    return f"{weekday}, {dt.day:02d} {month} {dt.year:04d} {dt.hour:02d}:{dt.minute:02d}:{dt.second:02d} GMT"


def count_decimal_digits(number):
    """
    Count the number of decimal digits in a number
    :param number: The number
    :return: The number of decimal digits
    """
    number_parts = repr(number).split('.')
    if len(number_parts) == 1:
        return 0
    else:
        if number_parts[1] == '0':
            return 0
        else:
            return len(number_parts[1].rstrip('0'))


def get_variable_bounds(m, data, basis):
    coeff_lb_count = defaultdict(int)
    coeff_ub_count = defaultdict(int)

    min_exponent = float("inf")
    max_exponent = -float("inf")

    for var in m.getVars():

        if var.LB != 0 and abs(var.LB) < gp.GRB.INFINITY:
            exponent = int(math.floor(math.log(abs(var.LB), basis)))
            coeff_lb_count[exponent] += 1
            if exponent < min_exponent: min_exponent = exponent
            if exponent > max_exponent: max_exponent = exponent

        if var.UB != 0 and abs(var.UB) < gp.GRB.INFINITY:
            exponent = int(math.floor(math.log(abs(var.UB), basis)))
            coeff_ub_count[exponent] += 1
            if exponent < min_exponent: min_exponent = exponent
            if exponent > max_exponent: max_exponent = exponent

    if min_exponent == float("inf"):
        return [], []

    else:

        coeff_lb_result = [[exponent, coeff_lb_count[exponent]] for exponent in range(min_exponent, max_exponent + 1)]
        coeff_ub_result = [[exponent, coeff_ub_count[exponent]] for exponent in range(min_exponent, max_exponent + 1)]
        return coeff_lb_result, coeff_ub_result


def get_obj_c_frequencies(model, basis):
    int_count = defaultdict(int)
    cont_count = defaultdict(int)
    int_result = []
    cont_result = []

    for x in model.getVars():
        if x.Obj != 0:
            if x.VType == 'C':
                cont_count[int(math.floor(math.log(abs(x.Obj), basis)))] += 1
            else:
                int_count[int(math.floor(math.log(abs(x.Obj), basis)))] += 1

    if len(cont_count) > 0:
        cont_result = [ [exponent, cont_count[exponent]] for exponent in range(min(cont_count.keys()), max(cont_count.keys()) + 1) ]

    if len(int_count) > 0:
        int_result = [ [exponent, int_count[exponent]] for exponent in range(min(int_count.keys()), max(int_count.keys()) + 1) ]

    return cont_result, int_result
