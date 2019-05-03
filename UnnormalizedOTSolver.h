#ifndef _UnnormalizedOTSolver_
#define _UnnormalizedOTSolver_

#include "double_array_1d.h"
#include "double_array_2d.h"
#include "double_array_3d.h"
#include<cassert>
#include<cmath>
#include<vector>

using DoubleArray1D = Math270A::DoubleArray1D;
using DoubleArray2D = Math270A::DoubleArray2D;
using DoubleArray3D = Math270A::DoubleArray3D;

class UnnormalizedOTSolver {
public:
  UnnormalizedOTSolver(size_t nt, size_t nx, size_t ny, double tau_1 = 1e-2, double tau_2 = 1e-1, double alpha = 1.0, bool enforce_zero_f = false);

  void Solve();

  void SetObjectiveValueStoppingCriterion(double val) { assert(val >= 0.0); objective_value_stopping_criterion = val; }
  void SetConstraintViolationValueStoppingCriterion(double val) { assert(val >= 0.0); constraint_violation_stopping_criterion = val; }
  void SetMaxIterationStoppingCriterion(unsigned int iter) { max_iteration_stopping_criterion = iter; }

  void UpdateValues();

  void UpdatePhi();
  void UpdateM();
  void UpdateU();
  void UpdateF();

  void NormalizeU();

  void SetInputs(const DoubleArray2D& mu_1, const DoubleArray2D& mu_2);

  double CheckConstraintViolation() const;

  double CheckObjectiveValue() const;
  double GetDerivStencilDualityGapTime() const;
  double GetGradDivStencilDualityGap();


  double GetDivMIntegral() const;

  /**************************************************/
  // Calculus
  /**************************************************/

  // Abstract
  double DtPhi(size_t k, size_t i, size_t j) const;
  double DxPhi(size_t k, size_t i, size_t j) const;
  double DyPhi(size_t k, size_t i, size_t j) const;

  double DtU(size_t k, size_t i, size_t j) const;
  double DtUTilde(size_t k, size_t i, size_t j) const;

  double DivM(size_t k, size_t i, size_t j) const;
  double DivMTilde(size_t k, size_t i, size_t j) const;

  double GetSpacialIntegral(const DoubleArray3D& ar, size_t k) const;


  // Translations

  double SpaceDerivPrimalX(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;
  double SpaceDerivDualX(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;

  double SpaceDerivPrimalY(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;
  double SpaceDerivDualY(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;

  double TimeDerivPrimal(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;
  double TimeDerivDual(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;


  // Numerical

  double TimeDerivCentralPrimal(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;
  double TimeDerivCentralDual(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;

  double TimeDerivForwardsPrimal(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;
  double TimeDerivForwardsDual(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;

  double SpaceDerivCentralPrimalX(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;
  double SpaceDerivCentralDualX(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;

  double SpaceDerivCentralPrimalY(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;
  double SpaceDerivCentralDualY(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const;

  double GetTCellWidth() const;
  double GetXCellWidth() const;
  double GetYCellWidth() const;

  DoubleArray3D m_x;
  DoubleArray3D m_y;
  DoubleArray3D phi;
  DoubleArray3D u;
  DoubleArray1D f;
  DoubleArray3D m_x_tilde;
  DoubleArray3D m_y_tilde;
  DoubleArray3D u_tilde;
  DoubleArray1D f_tilde;
  const size_t n_t;
  const size_t n_x;
  const size_t n_y;
  double tau_1;
  double tau_2;
  double alpha;

  bool enforce_zero_f;

  /*
    Iteration Parameters
  */
  // Stopping Criterion
  double objective_value_stopping_criterion = 0.0;
  double constraint_violation_stopping_criterion = 0.0;
  unsigned int max_iteration_stopping_criterion = UINT_MAX;

  // The largest value that norm(m_x) or norm(m_y) can obtain before stopping the iteration due to divergence.
  // A larger value means that the algorithm will take longer to report a divergence.
  double divergence_tolerance = 1e10;

  // The number of iterations that will happen minimum, even if a stopping criterion is met.
  unsigned int min_iterations = 2 * n_t;
  
  std::vector<double> objective_value;
  std::vector<double> constraint_violation;

  unsigned int current_iteration = 0;
  double current_objective_value = std::numeric_limits<double>::max();
  double current_constraint_violation= std::numeric_limits<double>::max();

};

class UnnormalizedOTSolverDivergenceError : public std::logic_error {
public:
  UnnormalizedOTSolverDivergenceError(const UnnormalizedOTSolver& solver);
  const char* what() const override;

  unsigned int iteration;
  double objective_value;
  double constraint_value;
};

std::ostream& operator<<(std::ostream& ostr, const UnnormalizedOTSolver& solver);

double GetRealRootsOfCubic(double a, double b, double c, double d);
#endif