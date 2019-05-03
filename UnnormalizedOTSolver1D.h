#ifndef _UnnormalizedOTSolver1D_
#define _UnnormalizedOTSolver1D_

#include "double_array_1d.h"
#include "double_array_2d.h"
#include<cassert>
#include<cmath>
#include<vector>

using DoubleArray1D = Math270A::DoubleArray1D;
using DoubleArray2D = Math270A::DoubleArray2D;

class UnnormalizedOTSolver1D {
public:
	UnnormalizedOTSolver1D(size_t nt, size_t nx, double tau_1 = 1e-2, double tau_2 = 1e-1, double alpha = 1.0, bool enforce_zero_f = false);

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

	void SetInputs(const DoubleArray1D& mu_1, const DoubleArray1D& mu_2);

	double CheckConstraintViolation() const;

	double CheckObjectiveValue() const;

  double GetDerivStencilDualityGapSpace();
  double GetDerivStencilDualityGapTime() const;

  double GetDivMIntegral() const;

	/**************************************************/
	// Calculus
	/**************************************************/

	// Abstract
	double DtPhi(size_t k, size_t i) const;
	double DxPhi(size_t k, size_t i) const;

	double DtU(size_t k, size_t i) const;
	double DtUTilde(size_t k, size_t i) const;

	double DivM(size_t k, size_t i) const;
	double DivMTilde(size_t k, size_t i) const;

	double GetSpacialIntegral(const DoubleArray2D& ar, size_t k) const;


	// Translations

  double SpaceDerivPrimal(const DoubleArray2D& ar, size_t k, size_t i) const;
  double SpaceDerivDual(const DoubleArray2D& ar, size_t k, size_t i) const;

  double TimeDerivPrimal(const DoubleArray2D& ar, size_t k, size_t i) const;
  double TimeDerivDual(const DoubleArray2D& ar, size_t k, size_t i) const;

	// Numerical

  double TimeDerivCentralPrimal(const DoubleArray2D& ar, size_t k, size_t i) const;
  double TimeDerivCentralDual(const DoubleArray2D& ar, size_t k, size_t i) const;
  
  double TimeDerivBackwardsPrimal(const DoubleArray2D& ar, size_t k, size_t i) const;
  double TimeDerivBackwardsDual(const DoubleArray2D& ar, size_t k, size_t i) const;

  double TimeDerivForwardsPrimal(const DoubleArray2D& ar, size_t k, size_t i) const;
  double TimeDerivForwardsDual(const DoubleArray2D& ar, size_t k, size_t i) const;
  
  double SpaceDerivBackwardsPrimal(const DoubleArray2D& ar, size_t k, size_t i) const;
  double SpaceDerivBackwardsDual(const DoubleArray2D& ar, size_t k, size_t i) const;

  double SpaceDerivForwardsPrimal(const DoubleArray2D& ar, size_t k, size_t i) const;
  double SpaceDerivForwardsDual(const DoubleArray2D& ar, size_t k, size_t i) const;

  double SpaceDerivCentralPrimal(const DoubleArray2D& ar, size_t k, size_t i) const;
  double SpaceDerivCentralDual(const DoubleArray2D& ar, size_t k, size_t i) const;

	double GetTCellWidth() const;
	double GetXCellWidth() const;

  DoubleArray2D phi_old;

	DoubleArray2D m_x;
	DoubleArray2D phi;
	DoubleArray2D u;
	DoubleArray1D f;
	DoubleArray2D m_x_tilde;
	DoubleArray2D u_tilde;
	DoubleArray1D f_tilde;
	const size_t n_t;
	const size_t n_x;
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
  double current_constraint_violation = std::numeric_limits<double>::max();
};

class UnnormalizedOTSolver1DDivergenceError : public std::logic_error {
public:
  UnnormalizedOTSolver1DDivergenceError(const UnnormalizedOTSolver1D& solver);
  const char* what() const override;

  unsigned int iteration;
  double objective_value;
  double constraint_value;
};

std::ostream& operator<<(std::ostream& ostr, const UnnormalizedOTSolver1D& solver);

double GetRealRootsOfCubic(double a, double b, double c, double d);

#endif