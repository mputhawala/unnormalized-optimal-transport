#include "UnnormalizedOTSolver1D.h"
#include "FunctionUtilities.h"
#include <math.h>
#include <random>
#include <algorithm>

UnnormalizedOTSolver1D::UnnormalizedOTSolver1D(size_t nt, size_t nx,
  double tau_1, double tau_2, double alpha, bool enforce_zero_f)
  : m_x(nt, nx + 1), phi(nt, nx), u(nt, nx), f(nt)
  , m_x_tilde(nt, nx + 1), u_tilde(nt, nx), f_tilde(nt)
  , n_t(nt), n_x(nx), tau_1(tau_1), tau_2(tau_2), alpha(alpha), enforce_zero_f(enforce_zero_f){
  assert(GetDerivStencilDualityGapSpace() < 1e-13*n_x);
  assert(GetDerivStencilDualityGapTime() < 1e-13*n_t);
}

/*
Attempt to solve for the minimal transport plan, subject to the given stopping criterion.
*/
void UnnormalizedOTSolver1D::Solve() {
  assert(objective_value_stopping_criterion > 1e-13 ||
    constraint_violation_stopping_criterion > 1e-13 ||
    max_iteration_stopping_criterion < UINT_MAX);
  while (current_iteration < min_iterations ||
    (current_objective_value > objective_value_stopping_criterion &&
      current_constraint_violation > constraint_violation_stopping_criterion &&
      current_iteration < max_iteration_stopping_criterion)) {
    UpdateValues();
  }
}

void UnnormalizedOTSolver1D::UpdateValues() {
  UpdateU();
  UpdateM();
  if (!enforce_zero_f) {
    UpdateF();
  }
  UpdatePhi();
  ++current_iteration;
  current_objective_value = CheckObjectiveValue();
  objective_value.push_back(current_objective_value);
  current_constraint_violation = CheckConstraintViolation();
  constraint_violation.push_back(current_constraint_violation);

  // If it looks like the algorithm is diverging, then throw an error
  if (m_x.NormInf() > divergence_tolerance) {
    throw UnnormalizedOTSolver1DDivergenceError(*this);
  }
}

void UnnormalizedOTSolver1D::UpdateM() {
  double u_average;
  for (size_t k = 0; k < n_t; ++k) {
    // On account of the zero flux conditions, m is forever zero at the boundary
    for (size_t i = 1; i < n_x; ++i) {
      u_average = (u(k, i) + u(k, i - 1)) / 2.0;
      m_x_tilde(k, i) = -m_x(k, i);
      m_x(k, i) = u_average / (u_average + tau_1) * (tau_1 * DxPhi(k, i) + m_x(k, i));
      m_x_tilde(k, i) += 2 * m_x(k, i);
    }
  }
}

void UnnormalizedOTSolver1D::UpdateU() {
  double a, b, c, d;
  double m_x_mean;
  for (size_t k = 1; k < n_t - 1; ++k) {
    for (size_t i = 0; i < n_x; ++i) {
      m_x_mean = (m_x(k, i + 1) + m_x(k, i)) / 2.0;
      a = 1;
      b = -(tau_1 * DtPhi(k, i) + u(k, i));
      c = 0;
      d = -tau_1 / 2.0 * (pow(m_x_mean, 2));

      u_tilde(k, i) = -u(k, i);

      u(k, i) = GetRealRootsOfCubic(a, b, c, d);
      u(k, i) = fmax(u(k, i), 0);

      u_tilde(k, i) += 2 * u(k, i);
    }
  }
}

void UnnormalizedOTSolver1D::UpdateF() {
  for (size_t k = 1; k < n_t - 1; ++k) {
    f_tilde(k) = -f(k);
    f(k) = alpha / (alpha + tau_1) *(tau_1 * GetSpacialIntegral(phi, k) + f(k));
    f_tilde(k) += 2 * f(k);
  }
}

void UnnormalizedOTSolver1D::UpdatePhi() {
  for (size_t k = 0; k < n_t; ++k) {
    for (size_t i = 0; i < n_x; ++i) {
      phi(k, i) += tau_2 * (DtUTilde(k, i) + DivMTilde(k, i) - f_tilde(k));
    }
  }
  /*for (size_t i = 0; i < n_x; ++i) {
    phi(0, i) += tau_2 * (DtUTilde(0, i) + DivMTilde(0, i) - f_tilde(0)) + (u(1, i) - u(0, i))*GetTCellWidth();
  }*/
}

void UnnormalizedOTSolver1D::NormalizeU() {
  double sum_rho_0 = 0;
  for (size_t i = 0; i < n_x; ++i) {
    sum_rho_0 += u(0, i);
  }
  for (size_t k = 1; k < n_t - 1; ++k) {
    double sum_u = 0, sum_u_tilde = 0;
    for (size_t i = 0; i < n_x; ++i) {
      sum_u += u(k, i);
      sum_u_tilde += u_tilde(k, i);
    }
    if (sum_u > 1e-1) {
      for (size_t i = 0; i < n_x; ++i) {
        u(k, i) *= sum_rho_0 / sum_u;
        u_tilde(k, i) *= sum_rho_0 / sum_u_tilde;
      }
    }
  }
}

double UnnormalizedOTSolver1D::CheckConstraintViolation() const {
  double violation_overall = 0;
  double violation_at_time;
  for (size_t k = 1; k < n_t - 1; ++k) {
    violation_at_time = 0;
    for (size_t i = 0; i < n_x; ++i) {
      violation_at_time += fabs(DtU(k, i) + DivM(k, i) - f(k));
    }
    violation_at_time *= GetXCellWidth() * GetTCellWidth();
    violation_overall += violation_at_time;
  }
  return violation_overall;
}

double UnnormalizedOTSolver1D::CheckObjectiveValue() const {
  double objective = 0;
  for (size_t k = 1; k < n_t - 1; ++k) {
    for (size_t i = 0; i < n_x; ++i) {
      if (fabs(pow(m_x(k, i), 2)) < 1e-8)
        continue;
      else
        objective += (pow(m_x(k, i), 2)) / (2 * u(k, i))
        * GetXCellWidth() * GetTCellWidth();
    }
    objective += pow(f(k), 2) * GetTCellWidth() / alpha;
  }
  return objective;
}

double UnnormalizedOTSolver1D::GetDerivStencilDualityGapSpace() {
  srand(n_x);
  for (size_t i = 0; i < n_x; ++i) {
    m_x(0, i) = double(rand()) / RAND_MAX;
    phi(0, i) = double(rand()) / RAND_MAX;
  }
  // zero flux condition
  m_x(0, 0) = 0;

  double integral_of_da_b = 0;
  for (size_t i = 1; i < n_x; ++i) {
    integral_of_da_b += DxPhi(0, i) * m_x(0, i);
  }
  integral_of_da_b *= GetXCellWidth();

  double integral_of_a_db = 0;
  for (size_t i = 0; i < n_x; ++i) {
    integral_of_a_db -= phi(0, i) * DivM(0, i);
  }
  integral_of_a_db *= GetXCellWidth();
  for (size_t i = 0; i < n_x; ++i) {
    m_x(0, i) = 0;
    phi(0, i) = 0;
  }
  return fabs(integral_of_a_db - integral_of_da_b);
}

double UnnormalizedOTSolver1D::GetDerivStencilDualityGapTime() const {
  DoubleArray2D a(n_t, 1), b(n_t, 1);
  srand(n_t);
  for (size_t k = 0; k < n_t; ++k) {
    a(k, 0) = double(rand()) / RAND_MAX;
    b(k, 0) = double(rand()) / RAND_MAX;
  }

  double integral_of_da_b = 0;
  for (size_t k = 0; k < n_t; ++k) {
    integral_of_da_b += TimeDerivPrimal(a, k, 0) * b(k, 0);
  }
  integral_of_da_b *= GetXCellWidth();

  double integral_of_a_db = 0;
  for (size_t k = 0; k < n_t; ++k) {
    integral_of_a_db -= a(k, 0) * TimeDerivDual(b, k, 0);
  }
  integral_of_a_db *= GetXCellWidth();
  return fabs(integral_of_a_db - integral_of_da_b);
}

double UnnormalizedOTSolver1D::GetDivMIntegral() const {
  DoubleArray2D div_m(n_t, n_x);
  for (size_t k = 0; k < n_t; ++k) {
    for (size_t i = 0; i < n_x; ++i) {
      div_m(k, i) = DivM(k, i);
    }
  }
  double max_so_far = 0;
  for (size_t k = 0; k < n_t; ++k) {
    max_so_far = fmax(max_so_far, GetSpacialIntegral(div_m, k));
  }
  return max_so_far;
}

/**************************************************/
// Calculus
/**************************************************/

// Abstract
inline double UnnormalizedOTSolver1D::DtPhi(size_t k, size_t i) const {
  return TimeDerivDual(phi, k, i);
}

inline double UnnormalizedOTSolver1D::DxPhi(size_t k, size_t i) const {
  return SpaceDerivDual(phi, k, i);
}

inline double UnnormalizedOTSolver1D::DivM(size_t k, size_t i) const {
  return SpaceDerivPrimal(m_x, k, i);
}

inline double UnnormalizedOTSolver1D::DivMTilde(size_t k, size_t i) const {
  return SpaceDerivPrimal(m_x_tilde, k, i);
}

inline double UnnormalizedOTSolver1D::DtU(size_t k, size_t i) const {
  return TimeDerivPrimal(u, k, i);
}

inline double UnnormalizedOTSolver1D::DtUTilde(size_t k, size_t i) const {
  return TimeDerivPrimal(u_tilde, k, i);
}

inline double UnnormalizedOTSolver1D::GetSpacialIntegral(const DoubleArray2D & ar, size_t k) const {
  double r = 0;
  for (size_t i = 0; i < n_x; ++i) {
    r += ar(k, i);
  }
  return r*GetXCellWidth();
}

// Translations

/*
Dirichlet (zero boundary) conditions on the primal variable.
*/
inline double UnnormalizedOTSolver1D::TimeDerivPrimal(const DoubleArray2D& ar, size_t k, size_t i) const {
  // return TimeDerivForwardsPrimal(ar, k, i);
  return TimeDerivCentralPrimal(ar, k, i);
  // return TimeDerivBackwardsPrimal(ar, k, i);
}

inline double UnnormalizedOTSolver1D::TimeDerivDual(const DoubleArray2D& ar, size_t k, size_t i) const {
  // return TimeDerivForwardsDual(ar, k, i);
  return TimeDerivCentralDual(ar, k, i);
  // return TimeDerivBackwardsDual(ar, k, i);
}

/*
Naumann (zero flux) conditions in the primal variable.
*/
inline double UnnormalizedOTSolver1D::SpaceDerivPrimal(const DoubleArray2D& ar, size_t k, size_t i) const {
  assert(i < n_x + 1);
  return (ar(k, i + 1) - ar(k, i)) / GetXCellWidth();
  // return SpaceDerivForwardsPrimal(ar, k, i);
  // return SpaceDerivCentralPrimal(ar, k, i);
  return SpaceDerivBackwardsPrimal(ar, k, i);
}

inline double UnnormalizedOTSolver1D::SpaceDerivDual(const DoubleArray2D& ar, size_t k, size_t i) const {
  assert(i > 0);
  return (ar(k, i) - ar(k, i - 1)) / GetXCellWidth();
  // return SpaceDerivForwardsDual(ar, k, i);
  // return SpaceDerivCentralDual(ar, k, i);
  return SpaceDerivBackwardsDual(ar, k, i);
}

// Numerical

inline double UnnormalizedOTSolver1D::TimeDerivCentralPrimal(const DoubleArray2D& ar, size_t k, size_t i) const {
  if (k > 0 && k < n_t - 1)
    return (ar(k + 1, i) - ar(k - 1, i)) / (2 * GetTCellWidth());
  else if (k == 0)
    return (ar(1, i) - ar(0, i)) / GetTCellWidth();
  else
    return (ar(n_t - 1, i) - ar(n_t - 2, i)) / GetTCellWidth();
}

inline double UnnormalizedOTSolver1D::TimeDerivCentralDual(const DoubleArray2D& ar, size_t k, size_t i) const {
  if (k > 1 && k < n_t - 2)
    return (ar(k + 1, i) - ar(k - 1, i)) / (2 * GetTCellWidth());
  else if (k == 0)
    return (ar(0, i) + ar(1, i) / 2) / GetTCellWidth();
  else if (k == 1)
    return (ar(2, i) / 2 - ar(0, i)) / GetTCellWidth();
  else if (k == n_t - 2)
    return (ar(n_t - 1, i) - ar(n_t - 3, i) / 2) / GetTCellWidth();
  else
    return (-ar(n_t - 1, i) - ar(n_t - 2, i) / 2) / GetTCellWidth();
}

inline double UnnormalizedOTSolver1D::TimeDerivBackwardsPrimal(const DoubleArray2D& ar, size_t k, size_t i) const {
  if (k == 0)
    return ar(0, i) / GetTCellWidth();
  else
    return (ar(k, i) - ar(k - 1, i)) / GetTCellWidth();
}

inline double UnnormalizedOTSolver1D::TimeDerivBackwardsDual(const DoubleArray2D& ar, size_t k, size_t i) const {
  if (k == n_t - 1)
    return -ar(n_t - 1, i) / GetTCellWidth();
  else
    return (ar(k + 1, i) - ar(k, i)) / GetTCellWidth();
}

inline double UnnormalizedOTSolver1D::TimeDerivForwardsPrimal(const DoubleArray2D& ar, size_t k, size_t i) const {
  if (k == n_t - 1)
    return (-ar(n_t - 1, i)) / GetTCellWidth();
  else
    return (ar(k + 1, i) - ar(k, i)) / GetTCellWidth();
}

inline double UnnormalizedOTSolver1D::TimeDerivForwardsDual(const DoubleArray2D& ar, size_t k, size_t i) const {
  if (k == 0)
    return ar(0, i) / GetTCellWidth();
  else
    return (ar(k, i) - ar(k - 1, i)) / GetTCellWidth();
}

inline double UnnormalizedOTSolver1D::SpaceDerivBackwardsPrimal(const DoubleArray2D& ar, size_t k, size_t i) const {
  return (ar(k, i) - ar(k, i - 1)) / GetXCellWidth();
  if (i == 0)
    return 0.0;
  else
    return (ar(k, i) - ar(k, i - 1)) / GetXCellWidth();
}

inline double UnnormalizedOTSolver1D::SpaceDerivBackwardsDual(const DoubleArray2D& ar, size_t k, size_t i) const {
  return (ar(k, i + 1) - ar(k, i)) / GetXCellWidth();
  if (i > 0 && i < n_x - 1)
    return (ar(k, i + 1) - ar(k, i)) / GetXCellWidth();
  else if (i == 0)
    return ar(k, 1) / GetXCellWidth();
  else
    return -ar(k, n_x - 1) / GetXCellWidth();
}

inline double UnnormalizedOTSolver1D::SpaceDerivForwardsPrimal(const DoubleArray2D& ar, size_t k, size_t i) const {
  if (i == n_x - 1)
    return 0.0;
  else
    return (ar(k, i + 1) - ar(k, i)) / GetXCellWidth();
}

inline double UnnormalizedOTSolver1D::SpaceDerivForwardsDual(const DoubleArray2D& ar, size_t k, size_t i) const {
  if (i > 0 && i < n_x - 1)
    return (ar(k, i) - ar(k, i - 1)) / GetXCellWidth();
  else if (i == 0)
    return ar(k, 0) / GetXCellWidth();
  else
    return -ar(k, n_x - 2) / GetXCellWidth();
}

inline double UnnormalizedOTSolver1D::SpaceDerivCentralPrimal(const DoubleArray2D& ar, size_t k, size_t i) const {
  if (i > 0 && i < n_x - 1)
    return (ar(k, i + 1) - ar(k, i - 1)) / (2 * GetXCellWidth());
  if (i == 0)
    return (ar(k, 1) - ar(k, 0)) / (2 * GetXCellWidth());
  else
    return (ar(k, n_x - 1) - ar(k, n_x - 2)) / (2 * GetXCellWidth());
}

inline double UnnormalizedOTSolver1D::SpaceDerivCentralDual(const DoubleArray2D& ar, size_t k, size_t i) const {
  if (i > 0 && i < n_x - 1)
    return (ar(k, i + 1) - ar(k, i - 1)) / (2 * GetXCellWidth());
  if (i == 0)
    return (ar(k, 1) + ar(k, 0)) / (2 * GetXCellWidth());
  else
    return (-ar(k, n_x - 1) - ar(k, n_x - 2)) / (2 * GetXCellWidth());
}

inline double UnnormalizedOTSolver1D::GetTCellWidth() const {
  return 1.0 / (n_t - 1.0);
}

inline double UnnormalizedOTSolver1D::GetXCellWidth() const {
  return 1.0 / (n_x - 1.0);
}

void UnnormalizedOTSolver1D::SetInputs(const DoubleArray1D& mu_1, const DoubleArray1D& mu_2) {
  assert(mu_1.GetIndex1Size() == mu_2.GetIndex1Size());

  m_x.SetToValue(0.0);
  phi.SetToValue(0.0);
  u.SetToValue(0.0);
  f.SetToValue(0.0);

  size_t last_ind = n_t - 1;
  for (size_t i = 0; i < mu_1.GetIndex1Size(); ++i) {
    u(0, i) = u_tilde(0, i) = mu_1(i);
    u(last_ind, i) = u_tilde(last_ind, i) = mu_2(i);
  }
}

std::ostream& operator<<(std::ostream& ostr, const UnnormalizedOTSolver1D& solver) {
  std::cout << "m_x = " << '\n' << solver.m_x << '\n';

  std::cout << "phi = " << '\n' << solver.phi << '\n';
  std::cout << "u = " << '\n' << solver.u << '\n';
  // std::cout << "u_tilde = " << '\n' << solver.u_tilde << '\n';
  std::cout << "f = " << '\n' << solver.f << '\n';

  std::cout << "tau_1 = " << solver.tau_1 << '\n';
  std::cout << "tau_2 = " << solver.tau_2 << '\n';
  // std::cout << "beta = " << solver.beta << '\n';
  return ostr;
}

UnnormalizedOTSolver1DDivergenceError::UnnormalizedOTSolver1DDivergenceError(const UnnormalizedOTSolver1D& solver)
  : logic_error("Iterative Algorithm has diverged.") {
  iteration = solver.current_iteration;
  objective_value = solver.current_objective_value;
  constraint_value = solver.current_constraint_violation;
}

const char* UnnormalizedOTSolver1DDivergenceError::what() const {
  return logic_error::what();
}