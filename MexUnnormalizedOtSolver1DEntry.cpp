#ifdef MATLAB_MEX_FILE
#include <iostream>
#include "double_array_1d.h"
#include "double_array_2d.h"
#include "double_array_3d.h"
#include "UnnormalizedOTSolver.h"
#include "UnnormalizedOTSolver1D.h"
#include "FunctionUtilities.h"

#if __cplusplus < 199711L
#define nullptr NULL
#endif

#include <math.h>
#include <matrix.h>
#include <mex.h>

using namespace std;

void SaveDoubleArray2DToMexArray(const DoubleArray2D& double_array_2d, mxArray*& mx_array_ptr) {
  size_t n_1 = double_array_2d.GetIndex1Size();
  size_t n_2 = double_array_2d.GetIndex2Size();
  int dims[3] = { static_cast<int>(n_1), static_cast<int>(n_2) };
  mx_array_ptr = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  auto data_ptr = mxGetPr(mx_array_ptr);
  for (size_t i_1 = 0; i_1 < n_1; ++i_1) {
    for (size_t i_2 = 0; i_2 < n_2; ++i_2) {
      data_ptr[i_1 + i_2*n_1] = double_array_2d(i_1, i_2);
    }
  }
}

void SaveDoubleArray1DToMexArray(const DoubleArray1D& double_array_1d, mxArray*& mx_array_ptr) {
  size_t n_1 = double_array_1d.GetIndex1Size();
  int dims[3] = { static_cast<int>(n_1) };
  mx_array_ptr = mxCreateNumericArray(1, dims, mxDOUBLE_CLASS, mxREAL);
  auto data_ptr = mxGetPr(mx_array_ptr);
  for (size_t i_1 = 0; i_1 < n_1; ++i_1) {
    data_ptr[i_1] = double_array_1d(i_1);
  }
}

void MexEchoInputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  mexPrintf("Echoing inputs\n");
  for (size_t i = 0; i < nrhs; i++) {
    const size_t n = mxGetDimensions(prhs[i])[0];
    const size_t m = mxGetDimensions(prhs[i])[1];
    double* arg = mxGetPr(prhs[i]);
    mexPrintf("Input %u:\n", i);
    mexPrintf("[");
    for (size_t i1 = 0; i1 < n; i1++) {
      mexPrintf("\t[");
      if (n > 1 && m > 1) {
        mexPrintf("(%u,%u)", i1, 0);
      }
      mexPrintf("%1.4f", arg[i1*m + 0]);
      if (m != 1) {
        mexPrintf(",\t");
      }
      for (size_t i2 = 1; i2 < m; i2++) {
        if (n > 1 && m > 1) {
          mexPrintf("(%u,%u)", i1, i2);
        }
        mexPrintf("%1.4f", arg[i1*m + i2]);
        if (i2 != m - 1) {
          mexPrintf(",\t");
        }
      }
      mexPrintf("]\n");
    }
    mexPrintf("]\n");
  }
}

void PrintUsage() {
  mexPrintf("Usage is: \n");
  mexPrintf("[double emd, 2darray m_x, 2darray phi, 2darray u, 1darray f, 2darray iteration & objective & constraint] = UnbalancedEMD(1darray rho_0, 1darray rho_1, int n_t, double tau_1, double tau_2, double alpha, bool enfore_zero_f, double objective_stop, double constraint_stop, int max_iterations);\n");
  return;
}
/*
Signature of the calling function should be:
[double emd, 2darray m_x, 2darray phi, 2darray u, 1darray f, 2darray iteration & objective & constraint] = UnbalancedEMD(1darray rho_0, 1darray rho_1, int n_t, double tau_1, double tau_2, double alpha, bool enfore_zero_f, double objective_stop, double constraint_stop, int max_iterations);
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
#ifdef _MATLAB_DEBUG
  MexEchoInputs(nlhs, plhs, nrhs, prhs);
  mexPrintf("Parsing inputs\n");
#endif
  if (nrhs != 10) {
    PrintUsage();
    mexErrMsgIdAndTxt("UnnormalizedOTSolver:InvalidArgument", "Wrong number of inputs");
  }
  if (nlhs != 6) {
    mexErrMsgIdAndTxt("UnnormalizedOTSolver:InvalidOutputs", "Wrong number of outputs");
  }
  for (size_t i = 0; i < 10; ++i) {
    assert(mxIsNumeric(phrs[i]));
    if (i < 2 && mxGetN(prhs[i]) != 1) {
      mexPrintf("Input %1u is not a vector\n", i);
      mexPrintf("Input %1u is: [%4u,%4u]\n", mxGetM(prhs[i]), mxGetN(prhs[i]));
      PrintUsage();
      mexErrMsgIdAndTxt("UnnormalizedOTSolver:InvalidArgument", "Inputs are the wrong size.");
    }
    if (i > 1 && (mxGetM(prhs[i]) != 1 || mxGetN(prhs[i]) != 1)) {
      mexPrintf("Input %1u is not scalar\n", i);
      mexPrintf("Input %1u is: [%4u,%4u]\n", mxGetM(prhs[i]), mxGetN(prhs[i]));
      PrintUsage();
      mexErrMsgIdAndTxt("UnnormalizedOTSolver:InvalidArgument", "Inputs are the wrong size.");
    }
  }

  const size_t n_x = static_cast<size_t>(mxGetDimensions(prhs[0])[0]);
  int      n_t = static_cast<int>(round(static_cast<double>(mxGetPr(prhs[2])[0])));
  double tau_1 = mxGetPr(prhs[3])[0];
  double tau_2 = mxGetPr(prhs[4])[0];
  double alpha = mxGetPr(prhs[5])[0];
  bool enforce_zero_f = (mxGetPr(prhs[6])[0] > 1e-12);
  double objective_value_stopping_val = mxGetPr(prhs[7])[0];
  double constraint_value_stopping_val = mxGetPr(prhs[8])[0];
  int max_iterations = static_cast<int>(mxGetPr(prhs[9])[0]);

  // Read rho0 and rho1
  DoubleArray1D rho0(n_x);
  DoubleArray1D rho1(n_x);
  double* rho0_pointer = mxGetPr(prhs[0]);
  double* rho1_pointer = mxGetPr(prhs[1]);

  for (size_t i = 0; i < n_x; ++i) {
    rho0(i) = rho0_pointer[i];
    rho1(i) = rho1_pointer[i];
  }

  // run solver
#ifdef _MATLAB_DEBUG
  mexPrintf("Capturing Parameters\n");
#endif
  UnnormalizedOTSolver1D solver(n_t, n_x, tau_1, tau_2, alpha, enforce_zero_f);
  solver.SetInputs(rho0, rho1);
  if (fabs(objective_value_stopping_val + 1) > 1e-13) {
#ifdef _MATLAB_DEBUG
    mexPrintf("Setting objective_stopping_val to %1.11f\n", objective_value_stopping_val);
#endif
    solver.SetObjectiveValueStoppingCriterion(objective_value_stopping_val);
  }
#ifdef _MATLAB_DEBUG
  else {
    mexPrintf("Not setting objective_stopping_val\n");
  }
#endif
  if (fabs(constraint_value_stopping_val + 1) > 1e-13) {
#ifdef _MATLAB_DEBUG
    mexPrintf("Setting constraint_value_stopping_val to %1.11f\n", constraint_value_stopping_val);
#endif
    solver.SetConstraintViolationValueStoppingCriterion(constraint_value_stopping_val);
  }
#ifdef _MATLAB_DEBUG
  else {
    mexPrintf("Not setting constraint_value_stopping_val\n");
  }
#endif
  if (fabs(max_iterations + 1) > 1e-13) {
#ifdef _MATLAB_DEBUG
    mexPrintf("Setting max_iterations to %9u\n", max_iterations);
#endif
    solver.SetMaxIterationStoppingCriterion(max_iterations);
  }
#ifdef _MATLAB_DEBUG
  else {
    mexPrintf("Not setting max_iterations\n");
  }
#endif
#ifdef _MATLAB_DEBUG
  mexPrintf("Applying Solver\n");
#endif
  try {
    solver.Solve();
  }
  catch (UnnormalizedOTSolverDivergenceError& er) {
    char message[100];
    sprintf(message, "Divergence on Iteration: %5u, Objective Val: %1.3f, Constraint Violation: %1.3f", er.iteration, er.objective_value, er.constraint_value);
    mexErrMsgIdAndTxt("UnnormalizedOTSolver:Divergence", message);
  }

  //for (int iteration = 1; iteration <= n_iters; ++iteration) {
  //  // cout << solver.CheckObjectiveValue() << "\t" << solver.CheckConstraintViolation() << endl;
  //  solver.UpdateValues();
  //}
  // set outputs
#ifdef _MATLAB_DEBUG
  mexPrintf("Setting Outputs\n");
  mexPrintf("Setting EMD\n");
#endif
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxGetPr(plhs[0])[0] = solver.CheckObjectiveValue();
#ifdef _MATLAB_DEBUG
  mexPrintf("Setting m_x\n");
#endif 
  SaveDoubleArray2DToMexArray(solver.m_x, plhs[1]);
#ifdef _MATLAB_DEBUG
  mexPrintf("Setting phi\n");
#endif
  SaveDoubleArray2DToMexArray(solver.phi, plhs[2]);
#ifdef _MATLAB_DEBUG
  mexPrintf("Setting u\n");
#endif
  SaveDoubleArray2DToMexArray(solver.u, plhs[3]);
#ifdef _MATLAB_DEBUG
  mexPrintf("Setting f\n");
#endif
  SaveDoubleArray1DToMexArray(solver.f, plhs[4]);
#ifdef _MATLAB_DEBUG
  mexPrintf("Setting objective and constraint\n");
#endif
  plhs[5] = mxCreateDoubleMatrix(solver.current_iteration, 3, mxREAL);
  double* optimization_output = mxGetPr(plhs[5]);
  const unsigned int n_iters = solver.current_iteration;
  for (unsigned int i = 0; i < n_iters; ++i) {
    optimization_output[i] = i;
  }
  for (unsigned int i = 0; i < n_iters; ++i) {
    optimization_output[n_iters + i] = solver.objective_value[i];
  }
  for (unsigned int i = 0; i < n_iters; ++i) {
    optimization_output[2 * n_iters + i] = solver.constraint_violation[i];
  }
#ifdef _MATLAB_DEBUG
  mexPrintf("Outputs Set.\n");
#endif
}

#endif