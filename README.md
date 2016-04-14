# SmoothedAggregationAlgebraicMultigrid.m
An aggregation-based algebraic multigrid matlab package. Includes sparsified smoothed aggregation.  

This is a Matlab software package for solving linear systems using aggregation based algebraic multigrid.  
Eran Treister and Irad Yavneh, 
Non-Galerkin Multigrid based on Sparsified Smoothed Aggregation. 
SIAM Journal on Scientific Computing, 37 (1), A30-A54, 2015.

The Technion---Israel Institute of Technology, Haifa 32000, Israel.
Contact email: eran@cs.technion.ac.il, irad@cs.technion.ac.il.

Please, cite the paper if you use our code.

Bug reports should be sent to eran@cs.technion.ac.il.

----------------------------------------------------------------------------

How to install this package:

Requirements: Matlab, compiler with OpenMP support.

1) Run makeAllMex.m in the Code/MEXfunc directory of this package.

2) Run demo.m and check that the example works.

----------------------------------------------------------------------------

How to use this package:

This package includes one main function SolveLinearSystem.m.

----------------------------------------------------------------------------
