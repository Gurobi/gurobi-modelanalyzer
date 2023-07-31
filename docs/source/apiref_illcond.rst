.. _APIRefLabel:

API Reference
#############


.. _APIkappa_explainLabel:

.. py:function:: model_analyzer.kappa_explain(model, data=None, KappaExact=-1, prmfile=None, relobjtype="LP", expltype="ROWS", method="DEFAULT", smalltol=1e-13, submatrix=False)

   Computes an explanation of ill conditioning for the basis matrix associated
   with the LP model in the first argument.   Writes the explanation to an LP
   or MPS file and returns a Gurobi model object to the subset of the basis
   matrix in the explanation.  Basis statuses associated with the model will
   be used without additional iterations; otherwise it will solve the LP to
   completion and use the resulting basis for the explanation computation.

   :param model: Required Gurobi model object containing an LP.
   :param data:  Do not provide.  Only for use by Gurobi routines.
   :param KappaExact: Optional for console display of condition number
		      information.  Set to 1 to obtain exact condition
		      number information (which can require computation time
		      on the order of 1/3 of the simplex iterations needed
		      to solve the LP).  Set to 0 to avoid this computation
		      and use a less computationally expensive estimate.  The
		      default of -1 allows the explainer to decide based on
		      problem size.
   :param prmfile:    Optional name of parameter file for the subproblem solved to derive the ill conditioning certificate.
   :param relobjtype: Optional type of subproblem to create for the calculating the ill conditioning certificate.  "LP" is the default, but "QP" specifies a quadratic objective.
   :param expltype:   Optional parameter to request an explanation by rows or columns.  The explanation consists of either a collection of rows or columns that are close to linearly dependent.  Specify "ROWS" (the default) or "COLS".  One type of explanation is often much smaller than the other.
   :param method:     Optional method used for explainer computation.  Non default settings may help generate smaller explanation.

		      "DEFAULT":   No regularization of subproblem objective

		      "ANGLES":    Perform simpler calculation based on inner products.  This can be potentially faster, but it only finds explanations of two rows or columns.

		      "LASSO":     One norm regularization of subproblem objective.

		      "RLS":       Two norm regularization of subproblem objective.
   :param smalltol:   Optional tolerance below which ill conditioning certificate values are all treated as zero, i.e., the row or column associated with the value is filtered out of the explanation.  If left at the default of 1e-13, row or column norms and machine precision will be incorporated into the filtering decision.
   :param submatrix:  Optional flag to try to postprocess the explanation down to a smaller submatrix.   Default is False.
   :return:           A Gurobi model object containing the basis matrix rows or columns in the explanation that was written to an LP or MPS file.

.. _APIangle_explainLabel:

.. py:function:: model_analyzer.angle_explain(model, howmany=1, partol=1e-6)

   Searches for ill conditioning explanations consisting of pairs of near
   parallel rows and near parallel columns of the basis matrix.  Basis statuses
   associated with the model will be used without additional iterations;
   otherwise it will solve the LP to completion and use the resulting basis
   for the explanation computation.

   :param model:      Required Gurobi model object containing an LP.
   :param howmany:    Optional number of near parallel row and near parallel columns to search for.
   :param partol:     Tolerance below which two row or column vectors are considered parallel.  Used as a relative tolerance when comparing the inner product of two vectors with the product of their L1 norms.  The default setting is 1e-6.
   :return:           A list of tuples of almost parallel rows and almost parallel columns, and the model associated with the basis matrix from which the list was derived.



.. _APImatrix_bitmapLabel:

.. py:function:: model_analyzer.matrix_bitmap(model)

   Utility routine to facilitate examination of explainer output by
   printing a bit map of the nonzero structure.   This can be used to
   examine the nonzero structure of the LP or MPS file containing the
   explanation, or it can be used to obtain bit maps for matrices of
   arbitrary models rather than the explainer model.  The matplotlib
   package must be installed; otherwise this function will return
   immediately.
  		 
   :param model:      Required Gurobi model object containing an LP.
   :return:           Nothing.  The bit map appears on the screen and
		      can be saved to a file.


.. _APIconverttofractionsLabel:

.. py:function:: model_analyzer.converttofractions(vals)

   Utility routine to convert decimals into the nearest rational
   approximation.  This is particularly useful for converting
   truncated or rounded repeating decimals that appear in an
   LP or MPS file generated by the
   :ref:`kappa_explain <APIkappa_explainLabel>` function.
  		 
   :param vals:       Required array of decimal values to convert.
   :return:           Nothing.  Results are printed to the screen.
		      
