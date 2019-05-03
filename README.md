# unnormalized-optimal-transport-
A repository containing code for computing theUunnormalized Wasserstein-2 (UW-2) distance in one and two dimensions.

The code is written in C++, with Mex wrappers so that the code can be called with Matlab. The files contents are as follows:

  UnnormalizedOTSolver1D.cpp/h - A C++ class that can be used to actually compute the UW-2 in 1 dimensions. For usage examples, see main.cpp

  UnnormalizedOTSolver.cpp/h - Same as above, but used for computing UW-2 in 2 dimensions.

  main.cpp - Includes both example usage for C++ interface as well as mex entry point for computing the UW-2 in 2 dimensions. If you compile the code without a mex compiler, the Matlab portion of the code will be automatically excluded. See mex_unbalanced_emd_code.m for mex example.
  
  MexUnnormalizedOtSolver1DEntry.cpp - Entry point for mex'ing code for computing the 1 dimensional UW-2. See mex_unbalanced_emd_code_1d.m for example usage.
  
  FunctionUtilities.cpp/h - Additional miscellaneous code.
  
  double_array_Nd.cpp/h - Data structues used to abstract a numerical grid in N dimensions. You can have the code use your favorite numerical grid class (provided that the interface is the same) by changing the using statement at the top of UnnormalizedOTSolver1D.h/UnnormalizedOTSolver.h.
  
  mex_unbalanced_emd_code.m - Example code that you can use to mex the UW_2 2 dimensional code so that you can call it from Matlab.
  
  mex_unbalanced_emd_code_1d.m - Same as above, but in 1 dimension
