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

#ifdef MATLAB_MEX_FILE
#include <math.h>
#include <matrix.h>
#include <mex.h>
#endif

using namespace std;
using Math270A::DoubleArray3D;

#ifdef MATLAB_MEX_FILE

void SaveDoubleArray3DToMexArray(const DoubleArray3D& double_array_3d, mxArray*& mx_array_ptr) {
  size_t n_1 = double_array_3d.GetIndex1Size();
  size_t n_2 = double_array_3d.GetIndex2Size();
  size_t n_3 = double_array_3d.GetIndex3Size();
  int dims[3] = { static_cast<int>(n_1), static_cast<int>(n_2), static_cast<int>(n_3) };
  mx_array_ptr = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  auto data_ptr = mxGetPr(mx_array_ptr);
  for (size_t i_1 = 0; i_1 < n_1; ++i_1) {
    for (size_t i_2 = 0; i_2 < n_2; ++i_2) {
      for (size_t i_3 = 0; i_3 < n_3; ++i_3) {
        data_ptr[i_1 + i_2*n_1 + i_3*n_1*n_2] = double_array_3d(i_1, i_2, i_3);
      }
    }
  }
}

void SaveDoubleArray2DToMexArray(const DoubleArray2D& double_array_2d, mxArray*& mx_array_ptr) {
  size_t n_1 = double_array_2d.GetIndex1Size();
  size_t n_2 = double_array_2d.GetIndex2Size();
  int dims[3] = { static_cast<int>(n_1), static_cast<int>(n_2)};
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
  int dims[3] = { static_cast<int>(n_1)};
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
  mexPrintf("[double emd, 3darray m_x, 3darray m_y, 3darray phi, 3darray u, 1darray f, 2darray iteration & objective & constraint] = UnbalancedEMD(2darray rho_0, 2darray rho_1, int n_t, double tau_1, double tau_2, double alpha, bool enfore_zero_f, double objective_stop, double constraint_stop, int max_iterations);\n");
  return;
}
/*
Signature of the calling function should be:
[double emd, 3darray m_x, 3darray m_y, 3darray phi, 3darray u, 1darray f, 2darray iteration & objective & constraint] = UnbalancedEMD(2darray rho_0, 2darray rho_1, int n_t, double tau_1, double tau_2, double alpha, bool enfore_zero_f, double objective_stop, double constraint_stop, int max_iterations);
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
#ifdef _MATLAB_DEBUG
  MexEchoInputs(nlhs, plhs, nrhs, prhs);
  mexPrintf("Parsing inputs\n");
#endif
  if (nrhs != 10) {
    PrintUsage();
    mexErrMsgIdAndTxt("UnnormalizedOTSolver:InvalidArgument","Wrong number of inputs");
  }
  if (nlhs != 7) {
    mexErrMsgIdAndTxt("UnnormalizedOTSolver:InvalidOutputs", "Wrong number of outputs");
  }
  for (size_t i = 0; i < 10; ++i) {
    assert(mxIsNumeric(phrs[i]));
    if (i > 1 && (mxGetM(prhs[i]) != 1 || mxGetN(prhs[i]) != 1)) {
      mexPrintf("Input %1u is not scalar\n", i);
      mexPrintf("Input %1u is: [%4u,%4u]\n", mxGetM(prhs[i]), mxGetN(prhs[i]));
      PrintUsage();
      mexErrMsgIdAndTxt("UnnormalizedOTSolver:InvalidArgument", "Inputs are the wrong size.");
    }
  }

  const size_t n_x = static_cast<size_t>(mxGetDimensions(prhs[0])[0]);
  const size_t n_y = static_cast<size_t>(mxGetDimensions(prhs[0])[1]);
  assert(mxGetDimensions(prhs[1])[0] == n_x);
  assert(mxGetDimensions(prhs[1])[1] == n_y);
  int      n_t = static_cast<int>(round(static_cast<double>(mxGetPr(prhs[2])[0])));
  double tau_1 = mxGetPr(prhs[3])[0];
  double tau_2 = mxGetPr(prhs[4])[0];
  double alpha = mxGetPr(prhs[5])[0];
  bool enforce_zero_f = (mxGetPr(prhs[6])[0] > 1e-12);
  double objective_value_stopping_val = mxGetPr(prhs[7])[0];
  double constraint_value_stopping_val = mxGetPr(prhs[8])[0];
  int max_iterations = static_cast<int>(mxGetPr(prhs[9])[0]);

  // Read rho0 and rho1
  DoubleArray2D rho0(n_x, n_y);
  DoubleArray2D rho1(n_x, n_y);
  double* rho0_pointer = mxGetPr(prhs[0]);
  double* rho1_pointer = mxGetPr(prhs[1]);

  for (size_t i = 0; i < n_x; ++i) {
    for (size_t j = 0; j < n_y; ++j) {
      // Matlab transposes inputs
      rho0(i, j) = rho0_pointer[i*n_y + j];
      rho1(i, j) = rho1_pointer[i*n_y + j];
    }
  }

  // run solver
#ifdef _MATLAB_DEBUG
  mexPrintf("Capturing Parameters\n");
#endif
  UnnormalizedOTSolver solver(n_t, n_x, n_y);
  solver.SetInputs(rho0, rho1);
  solver.tau_1 = tau_1;
  solver.tau_2 = tau_2;
  solver.alpha = alpha;
  solver.enforce_zero_f = enforce_zero_f;
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
  SaveDoubleArray3DToMexArray(solver.m_x, plhs[1]);
#ifdef _MATLAB_DEBUG
  mexPrintf("Setting m_y\n");
#endif
  SaveDoubleArray3DToMexArray(solver.m_y, plhs[2]);
#ifdef _MATLAB_DEBUG
  mexPrintf("Setting phi\n");
#endif
  SaveDoubleArray3DToMexArray(solver.phi, plhs[3]);
#ifdef _MATLAB_DEBUG
  mexPrintf("Setting u\n");
#endif
  SaveDoubleArray3DToMexArray(solver.u, plhs[4]);
#ifdef _MATLAB_DEBUG
  mexPrintf("Setting f\n");
#endif
  SaveDoubleArray1DToMexArray(solver.f, plhs[5]);
#ifdef _MATLAB_DEBUG
  mexPrintf("Setting objective and constraint\n");
#endif
  plhs[6] = mxCreateDoubleMatrix(solver.current_iteration, 3, mxREAL);
  double* optimization_output = mxGetPr(plhs[6]);
  const unsigned int n_iters = solver.current_iteration;
  for (unsigned int i = 0; i < n_iters; ++i) {
    optimization_output[i] = i;
  }
  for (unsigned int i = 0; i < n_iters; ++i) {
    optimization_output[n_iters + i] = solver.objective_value[i];
  }
  for (unsigned int i = 0; i < n_iters; ++i) {
    optimization_output[2*n_iters + i] = solver.constraint_violation[i];
  }
#ifdef _MATLAB_DEBUG
  mexPrintf("Outputs Set.\n");
#endif
}

#else

double Index2Coordinate(const DoubleArray1D& ar, size_t index) {
  return double(index) / (ar.GetIndex1Size() - 1);
}

double Index2Coordinate(const DoubleArray2D& ar, size_t index, size_t dim = 1) {
  if (dim == 1)
    return double(index) / (ar.GetIndex1Size() - 1);
  else
    return double(index) / (ar.GetIndex2Size() - 1);
}

size_t Coordinate2Index(const DoubleArray1D& ar, double coordinate) {
  return size_t (round(abs(coordinate * (ar.GetIndex1Size() - 1))));
}

void Initialize1DDoubleArrayToGaussian(double mu, double sigma, DoubleArray1D& ar) {
  double x;
  double total_mass = 0;
  for (size_t i = 0; i < ar.GetIndex1Size(); ++i) {
    x = Index2Coordinate(ar, i);
    ar(i) = exp(-pow(x - mu, 2) / (2.0 * pow(sigma, 2)));
    total_mass += ar(i);
  }
  ar /= total_mass / (ar.GetIndex1Size() - 1);
}

void Initialize2DDoubleArrayToGaussian(double mu_1, double mu_2, double sigma_1, double sigma_2, DoubleArray2D& ar) {
  double x, y;
  double total_mass = 0;
  for (size_t i = 0; i < ar.GetIndex1Size(); ++i) {
    for (size_t j = 0; j < ar.GetIndex2Size(); ++j) {
      x = Index2Coordinate(ar, i, 1);
      y = Index2Coordinate(ar, j, 2);
      ar(i, j) = exp(-pow(x - mu_1, 2) / (2.0 * pow(sigma_1, 2)) - pow(y - mu_2, 2) / (2.0 * pow(sigma_2, 2)));
      total_mass += ar(i, j);
    }
  }
  ar /= total_mass /( (ar.GetIndex1Size() - 1) * (ar.GetIndex2Size() - 1));
}

void InitDoubleArray(DoubleArray3D& in) {
  size_t val = 0;
  for (size_t i = 0; i < in.GetIndex1Size(); ++i) {
    for (size_t j = 0; j < in.GetIndex2Size(); ++j) {
      for (size_t k = 0; k < in.GetIndex3Size(); ++k) {
        in(i, j, k) = double(val++);
      }
    }
  }
}

void CreateGaussianDensities(DoubleArray2D& rho0, DoubleArray2D& rho1) {
  size_t n1 = rho0.GetIndex1Size();
  size_t n2 = rho0.GetIndex2Size();
  double sum0 = 0, sum1 = 0;

  for (size_t i = 0; i < n2; i++) {
    for (size_t j = 0; j < n1; j++) {
      double x = j / (n1*1.0);
      double y = i / (n2*1.0);

      double r1 = sqrt((x - 1.0 / 3)*(x - 1.0 / 3) + (y - 1.0 / 3)*(y - 1.0 / 3));
      double r2 = sqrt((x - 2.0 / 3)*(x - 2.0 / 3) + (y - 2.0 / 3)*(y - 2.0 / 3));

      rho0(i, j) = exp(-r2*r2 * 2);
      rho1(i, j) = exp(-r1*r1 * 2);
      sum0 += rho0(i, j);
      sum1 += rho1(i, j);


    }
  }
  sum0 /= n1*n2;
  sum1 /= n1*n2;
  for (size_t i = 0; i < n1; i++) {
    for (size_t j = 0; j < n2; ++j) {
      rho0(i, j) /= sum0;
      rho1(i, j) /= sum1;
    }
  }
}

void RunDoubleArray3DUnitTests() {
  // empty array construction
  DoubleArray3D array;
  cout << "Empty Array: " << array << endl;

  array.Resize(1, 1, 1);
  cout << "1-1-1 Array: " << array << endl;

  array.Resize(2, 2, 2);
  cout << "2-2-2 Array: " << array << endl;

  cout << "array.GetIndex1Size() = " << array.GetIndex1Size() << endl;
  cout << "array.GetIndex2Size() = " << array.GetIndex2Size() << endl;
  cout << "array.GetIndex3Size() = " << array.GetIndex3Size() << endl;


  // Comumn major labeling
  InitDoubleArray(array);
  cout << "2-2-2 Array after initialization: " << array << endl;

  // Arithmetic
  array *= 2;
  cout << "Array after array *= 2: " << array << endl;

  array /= 2;
  cout << "Array after array *= 2: " << array << endl;

  cout << "Array Inf Norm = " << array.NormInf() << endl;

  DoubleArray3D array2(array);
  cout << "Copy of Array: " << array2 << endl;
  array2 += array;
  array *= 2;
  array2 == array;
  cout << "Array + Array == 2*Array:" << (array2 == array) << endl;
}

void Run2DSolverTest() {
  size_t n_x = 35;
  size_t n_y = 35;
  size_t n_t = 15;
  double tau_1 = 1e-4;
  double tau_2 = 1e-1;
  double alpha = 1.0;
  int iters = 20000;
  UnnormalizedOTSolver solver(n_t, n_x, n_y, tau_1, tau_2, alpha, true);
  //output solver before anything is done
  //cout << solver << endl;

  // Set up the problem
  DoubleArray2D mu_1(n_x, n_y);
  DoubleArray2D mu_2(n_x, n_y);

  //CreateGaussianDensities(mu_1, mu_2);
  //for (size_t i = 0; i < n_x / 4; ++i) {
  //  for (size_t j = 0; j < n_y / 4; ++j) {
  //    mu_1(0, 0) = 1;
  //    mu_2(n_x - 1, n_y - 1) = 1;
  //    //mu_1(i, j) = 1; 16.0 / (n_x * n_y);
  //    //mu_2(n_x - 1 - i, n_y - 1 - j) = 1; 16.0 / (n_x * n_y);
  //  }
  //}
  mu_1(0, 0) = mu_1(1, 1) = mu_1(1, 0) = mu_1(0, 1) = 1;
  mu_2(n_x - 1, n_y - 1) = mu_2(n_x - 2, n_y - 2) = mu_2(n_x - 2, n_y - 1) = mu_2(n_x - 1, n_y - 2) = 1;
  solver.SetInputs(mu_1, mu_2);

  cout << "Mu 1:" << endl;
  cout << mu_1 << endl;

  cout << "Mu 2:" << endl;
  cout << mu_2 << endl;
  
  // If we know the exact solution in terms of u
  for (size_t i = 0; i < n_t; ++i)
    solver.u(i, i, i) = solver.u_tilde(i, i, i) = 
    solver.u(i, i, i + 1) = solver.u_tilde(i, i, i + 1) =
    solver.u(i, i + 1, i) = solver.u_tilde(i, i + 1, i) = 
    solver.u(i, i + 1, i + 1) = solver.u_tilde(i, i + 1, i + 1) = 1;

  cout << "Before Any Iterations:" << endl;
  cout << solver << endl;

  cout << "Objective/Constaint Value" << endl;
  for (int i = 1; i <= iters; ++i) {
    solver.UpdateValues();
    cout << solver.CheckObjectiveValue() << "\t" << solver.CheckConstraintViolation() << endl;
    //cout << "After " << i << " iterations" << endl;
    //cout << solver << endl;
  }
  cout << "Final Iteration:" << endl;
  cout << "Objective/Constaint Value" << endl;
  cout << solver.CheckObjectiveValue() << "\t" << solver.CheckConstraintViolation() << endl;
  cout << solver << endl;
}

void Run1DSolverTests() {
	size_t n_x = 7;
	size_t n_t = 5;
	int iter_test = 5000;
	DoubleArray1D rho1(n_x), rho2(n_x);

  //for (size_t i = 0; i < (n_x + 1) / 2; ++i) {
  //  rho1(i) = 1;
  //  rho2(n_x - i - 1) = 1;
  //}
  rho1(0) = rho1(1) = rho1(2) = 1;
  rho2(n_x - 1) = rho2(n_x - 2) = rho2(n_x - 3) = 1;

	UnnormalizedOTSolver1D solver(n_t, n_x);
  cout << "Space Derivative Duality Gap: " << solver.GetDerivStencilDualityGapSpace() << endl;
  cout << "Time Derivative Duality Gap: " << solver.GetDerivStencilDualityGapTime() << endl;
  
	solver.SetInputs(rho1, rho2);

  for (size_t i = 0; i < n_t; ++i) {
    // solver.u(i, i) = solver.u_tilde(i, i) = solver.u(i, i + 1) = solver.u_tilde(i, i + 1) = 1;
  }

	// cout << "Before solving anything" << endl;
	// cout << solver << endl;

	//cout << "Objective/Constaint Value/Integral of Div of M" << endl;
	 cout << "Objective/Constaint Value" << endl;
	for (int i = 1; i <= iter_test; ++i) {
		cout << solver.CheckObjectiveValue() << "\t" << solver.CheckConstraintViolation() << endl;
		solver.UpdateValues();
		// cout << solver << endl;
	}
  cout << "After " << iter_test << " iterations." << endl;
	cout << solver << endl;
}

void Run1DGaussianTests() {
  size_t n_t = 10;
  size_t n_x = 35;
  double tau_1 = 1e-3;
  double tau_2 = 1e-1;
  double alpha = 1.0;
  int iters = 20000;
  DoubleArray1D rho0(n_x), rho1(n_x);

  Initialize1DDoubleArrayToGaussian(0.0, .100, rho0);
  Initialize1DDoubleArrayToGaussian(2.0/3.0, .100, rho1);

  UnnormalizedOTSolver1D solver(n_t, n_x);
  solver.SetInputs(rho0, rho1);
  solver.tau_1 = tau_1;
  solver.tau_2 = tau_2;
  solver.alpha = alpha;

  cout << solver << endl;

  for (int i = 1; i <= iters; ++i) {
    cout << solver.CheckObjectiveValue() << "\t" << solver.CheckConstraintViolation() << endl;
    solver.UpdateValues();
  }
  cout << "After " << iters << " iterations." << endl;
  cout << solver << endl;
  cout << "Div M: " << solver.GetDivMIntegral() << endl;
}

void Run2DGaussianTests() {
  size_t n_t = 5;
  size_t n_x = 15;
  size_t n_y = 15;
  double tau_1 = 5e-3;
  double tau_2 = 1e-1;
  double alpha = 1.0;
  int iters = 200;
  bool force_f_to_zero = false;
  DoubleArray2D rho0(n_x, n_y), rho1(n_x, n_y);

  Initialize2DDoubleArrayToGaussian(0.0, 0.0, 0.10, 0.10, rho0);
  Initialize2DDoubleArrayToGaussian(1.0, 1.0, 0.10, 0.10, rho1);

  UnnormalizedOTSolver solver(n_t, n_x, n_y, tau_1, tau_2, alpha, force_f_to_zero);
  solver.SetInputs(rho0, rho1);
  solver.SetMaxIterationStoppingCriterion(iters);
  solver.Solve();
  // cout << solver << endl;
  /*
  for (int i = 1; i <= iters; ++i) {
    cout << solver.CheckObjectiveValue() << "\t" << solver.CheckConstraintViolation() << endl;
    solver.UpdateValues();
  }
  */
  cout << "After " << solver.current_iteration << " iterations." << endl;
  cout << solver << endl;
  cout << "Objective: " << solver.CheckObjectiveValue() << endl;
}

void RunCalculusTests() {
  const size_t n_x = 3;
  const size_t n_y = 3;
  UnnormalizedOTSolver solver(3, n_x, n_y);
  DoubleArray2D phi(n_x, n_y);
  for (size_t i = 0; i < n_x; ++i) {
    for (size_t j = 0; j < n_y; ++j) {
      phi(i, j) = 1.0 / (i + j + 1);
    }
  }

  for (size_t i = 0; i < n_x; ++i) {
    for (size_t j = 0; j < n_y; ++j) {
      solver.phi(1, i, j) = phi(i, j);
    }
  }
  cout << "Phi:" << endl;
  cout << phi;

  DoubleArray2D grad_phi_x(n_x, n_y), grad_phi_y(n_x, n_y);
  for (size_t i = 0; i < n_x; ++i) {
    for (size_t j = 0; j < n_y; ++j) {
      grad_phi_x(i, j) = solver.DxPhi(1, i, j);
      grad_phi_y(i, j) = solver.DyPhi(1, i, j);
    }
  }
  cout << "DxPhi:\n";
  cout << grad_phi_x;
  cout << "\nDyPhi:\n";
  cout << grad_phi_y;

  DoubleArray2D div_grad_phi(n_x, n_y);
  for (size_t i = 0; i < n_x; ++i) {
    for (size_t j = 0; j < n_y; ++j) {
      solver.m_x(1, i, j) = grad_phi_x(i, j);
      solver.m_y(1, i, j) = grad_phi_y(i, j);
    }
  }
  for (size_t i = 0; i < n_x; ++i) {
    for (size_t j = 0; j < n_y; ++j) {
      div_grad_phi(i, j) = solver.DivM(1, i, j);
    }
  }
  cout << "Div Grad Phi:\n";
  cout << div_grad_phi;
}

void CubicSovlerTest() {
  double a, b, c, d;
  a = 1;
  b = -2;
  c = 0;
  d = -6;
  cout << GetRealRootsOfCubic(a, b, c, d) << endl;
}

void Run1DTwoGaussiansToOne() {
  size_t n_t = 10;
  size_t n_x = 35;
  double tau_1 = 1e-3;
  double tau_2 = 1e-1;
  double alpha = 1.0;
  int iters = 2000;
  bool force_f_to_zero = false;
  DoubleArray1D rho_0(n_x), rho_1(n_x);

  Initialize1DDoubleArrayToGaussian(0.0, 0.1, rho_0);
  Initialize1DDoubleArrayToGaussian(0.33, 0.1, rho_1);
  rho_0 += rho_1;
  Initialize1DDoubleArrayToGaussian(1.0, 0.1, rho_1);

  UnnormalizedOTSolver1D solver(n_t, n_x, tau_1, tau_2, alpha, force_f_to_zero);
  solver.SetInputs(rho_0, rho_1);
  solver.SetMaxIterationStoppingCriterion(iters);

  solver.Solve();

  cout << solver << endl;
}


int main(int argc, char * argv[]) {
  // cout << "Hello World!" << endl;
  // RunDoubleArray3DUnitTests();
  // CubicSovlerTest();
  // RunCalculusTests();
  Run2DSolverTest();
  // Run1DSolverTests();
  // Run2DGaussianTests();
  // Run1DGaussianTests();
  //Run1DTwoGaussiansToOne();

  return 0;
}
#endif