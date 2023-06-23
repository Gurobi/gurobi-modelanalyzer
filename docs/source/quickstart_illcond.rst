.. _QSGuideLabel:

Quick Start Guide
#################

.. _QSIntroductionLabel:

Introduction
************

The Ill Conditioning Explainer can be used to obtain an explanation
of ill conditioned basis matrices for linear programs (LP).  It can be used
on the initial relaxation of a mixed integer program (MIP) by invoking the
Model.relax() function on the MIP to obtain the LP relaxation.  Other
problem types are not currently supported, although models with convex
quadratic objective functions can be turned into related LPs by replacing
the convex quadratic objective expression with a linear expression (including
just removing the quadratic part).  The currently implemented approach
uses the fact that an ill conditioned square matrix has linear combination
of a subset of rows or columns that is close to 0.  Mathematically,
for a square matrix B, itcalculates a nonzero linear combination vector y such
that :math:`||B^{T}y||_{\infty} \leq \epsilon` and :math:`||y|| >> \epsilon`.
The vector y can be viewed as a certificate of ill conditioning.  The elements
of y that are zero correspond to rows of the matrix B that do not
contribute to the ill conditioning.   The Ill Conditioning Explainer filters
out the rows corresponding to the zero elements of y.   The remaining rows
provide an explanation of the ill conditioning.   One can similarly calculate
a nonzero linear combination of the columns of B and derive subset of the
columns of B to explain the ill conditioning.   The Ill Conditioning
Explainer can do both.   For the special case of ill conditioning due to
pairs of almost parallel rows or columns, a separate method is available.


Running the Explainer
*********************

This section describes the 3 most basic ways to run the explainer.
Additional methods will be described later in the :doc:`advanced_usage_illcond`
section. The :ref:`kappa_explain <APIkappa_explainLabel>` method looks for
explanations of
arbitrary size. By default it will provide a row-based explanation, but a
column-based explanation can be obtained by setting the expltype function
argument to "COLS".   Here is the simplest sequence of commands to do this.


.. code-block:: python
		
   import gurobipy as gp
   import model_analyzer as ma

   m=gp.read("myillconditionedmodel.mps")
   m.optimize()
   kappa_explain(m)             # row-based; could also specify expltype="ROWS"


If a column-based explanation is preferred, call

.. code-block:: python

   kappa_explain(m, expltype="COLS")

One can request both types of explanations by consecutive kappa_explain calls
in a single code fragment.

The angle_explain method can be called with just the model argument, in which
case it will stop and return the first pair of almost parallel rows or columns.
Or it can be called with a numerical argument specifying the number of pairs
of near parallel rows and columns to return.


Interpreting the Output
***********************

The kappa_explain method outputs a file of the rows or columns of the basis
matrix that explain the ill conditioning using the most suitable problem
file format.  Row based explanations are written in LP format, while column
based explanations are written in MPS format.   The file name consists of
the name of the model as specified in the gurobipy Model.modelName attribute,
suffixed by "_kappaexplain", followed by the file format (.lp or .mps).
The method also returns the model associated with the LP or MPS file that
was written.

The LP format explanation consists first of a row consisting of the
computed linear combination of the basis matrix rows that takes on
small values.  Subsequent rows correspond to the constraints that are
involved in the ill conditioning.  The row names in the LP file
consist of the constraint name prefixed by a text string that
specifies the value of the element in the y vector that comprises the
certificate of infeasibility mentioned in the :ref:`QSIntroductionLabel`
section.  These are the row multipliers of the linear combination that
results in the first row of the LP file, which has the name
"GRB_Combined_Row".  The rows are listed in descending order by
largest absolute multiplier.   Here is an example of the constraints in
such an LP file::

  GRB_Combined_Row: = -1.75613e-08
  (mult=0.5000000001354858)R09bad: - 0.999999999 X02 + X01 - X03 = 0
  (mult=0.49999999963548586)R09: X02 - X01 + X03 = 0
  (mult=-4.999999859945595e-10)X46: - X03 + 0.109 X22 <= 0
  (mult=4.716981000007684e-10)R10: - 1.06 X01 + X04 = 0
  (mult=4.716981000007684e-10)R20: - 0.43 X22 + X26 = 0
  (mult=-4.716981000007684e-10)X50: X04 + X26 <= 310
  (mult=2.573301814737374e-10)X27: X22 <= 500


The last 5 rows all consist of small multipliers and constraint
coefficients on the order of 1 or less, so their contribution to the
explanation is minimal.  The first two basis rows provide the
explanation, and examination of these reveals that they intersect the
same 3 variables, and that the rows are almost parallel.  This is thie
cause of the ill conditioning.  Remedies depend on the model in
general, but in this case the modeler should focus on whether the
coefficient of 0.999999999 in constraint R09bad has a meaningful
difference fromn the coefficient in constraint R09, or if the
difference involves roundoff error in the data computation for the
model.  In the former case, the model needs to reflect that difference
in a way that it is larger; in the latter case the coefficient should
be cleaned up to the true value of 1.0.

The column based explanation is analagous, with an MPS file consisting
of a combined column named "GRB_Combined_Column" followed by the individual
columns in the explanation.   Here is the explanation for the same model::

  COLUMNS
    GRB_Combined_Column  R09       1.1705685309948421e-09
    (mult=1.240802673407655)X04  R10       1
    (mult=1.240802673407655)X04  X50       1
    (mult=-1.240802673407655)GRBslack_X50  X50       1
    (mult=1.1705685609891108)X02  R09       1
    (mult=1.1705685609891108)X02  X21       -1
    (mult=1.1705685609891108)X02  R09bad    -0.999999999
    (mult=1.1705685598185422)X01  R09       -1
    (mult=1.1705685598185422)X01  R10       -1.06
    (mult=1.1705685598185422)X01  X05       1
    (mult=1.1705685598185422)X01  R09bad    1
    (mult=1.1705685598185422)X01  X48       0.301
    (mult=-1.1705685598185422)GRBslack_X05  X05       1
    (mult=0.8862876247488983)X16  R13       1
    (mult=0.8862876247488983)X16  X51       1
    (mult=-0.8862876247488983)GRBslack_X51  X51       1
    (mult=0.8361204007065078)X14  X21       1.4
    (mult=0.8361204007065078)X14  R12       1
    (mult=0.8361204007065078)X06  R12       -1
    (mult=0.8361204007065078)X06  R13       -1.06
    (mult=0.8361204007065078)X06  X17       1
    (mult=0.8361204007065078)X06  X49       0.301
    (mult=-0.8361204007065078)GRBslack_X17  X17       1
    (mult=-0.5033444809736455)GRBslack_X49  X49       1
    (mult=-0.35234113650538124)X23  R19       1
    (mult=-0.35234113650538124)X23  X44       -1
    (mult=0.35234113650538124)X24  R19       1
    (mult=0.35234113650538124)X24  X48       -1
    (mult=-0.2516722403609866)X36  X44       1.4
    (mult=-0.2516722403609866)X36  R23       -1
    (mult=-0.2516722403609866)X37  R23       1
    (mult=-0.2516722403609866)X37  X49       -1
  RHS

For this model, the row-based explanation is easier to interpret, so
we will not examine the column-based explanation in detail.  However,
note that the problematic coefficient of 0.999999999 in constraint
R09bad does appear in the column-based output.

The angle_explain method does not output an LP or MPS file containing
the basis matrix rows or columns that explain the ill conditioning,
as it is capable of providing multiple (simple) explanations at once.
Rather, it returns separate tuplelists of pairs of almost parallel rows and
almost parallel columns, followed by the model associated with the basis
matrix that generated those tuplelists.   Note that one or both tuplelists
may be empty if no almost parallel rows or columns were detected.
Running the angle_explain method on the same model that generated the
previous LP and MPS files::

   >>> ma.angle_explain(m)
       ([(<gurobi.Constr R09>, <gurobi.Constr R09bad>)], [],
       <gurobi.Model Continuous instance basismodel: 28 constrs, 28 vars,
       No parameter changes>)
   >>> 

Thus, it detected the same two constraints as the row-based explanation,
and found no almost parallel columns.


Suggested Usage Quick Start
***************************

Explainer output may be small and straightforward to interpret, or it may
consists of hundreds or thousands of constraints or variables when run on
large models.   A model detailed discussion will appear in the
:doc:`advanced_usage_illcond` section regarding how to interpret large
explanations. However, the recommended approach when getting started
is to first request
a row-based explanation, but if it is too large to interpret then request
the column based explanation.   In many cases, one of the two explanations
may be much smaller and easier to interpret than the other.   If both are
large, try the angle_explain routine.   If none of these approach yield
anything, look at the :doc:`advanced_usage_illcond` section for additional
information on how to interpret the output.
