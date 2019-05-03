#include "UnnormalizedOTSolver.h"
#include "FunctionUtilities.h"

UnnormalizedOTSolver::UnnormalizedOTSolver(size_t nt, size_t nx, size_t ny,
  double tau_1, double tau_2, double alpha, bool enforce_zero_f)
  : m_x(nt, nx + 1, ny), m_y(nt, nx, ny + 1), phi(nt, nx, ny), u(nt, nx, ny), f(nt)
  , m_x_tilde(nt, nx + 1, ny), m_y_tilde(nt, nx, ny + 1), u_tilde(nt, nx, ny), f_tilde(nt)
  , n_t(nt), n_x(nx), n_y(ny), tau_1(tau_1), tau_2(tau_2), alpha(alpha), enforce_zero_f(enforce_zero_f){
  assert(GetDerivStencilDualityGapTime() < 1e-13*n_t);
  assert(GetGradDivStencilDualityGap() < 1e-13*n_x*n_y);
}

/*
  Attempt to solve for the minimal transport plan, subject to the given stopping criterion.
*/
void UnnormalizedOTSolver::Solve() {
  assert(objective_value_stopping_criterion > 1e-13 ||
    constraint_violation_stopping_criterion > 1e-13 ||
  max_iteration_stopping_criterion < UINT_MAX);
  while ( current_iteration < min_iterations ||
    (current_objective_value > objective_value_stopping_criterion &&
    current_constraint_violation > constraint_violation_stopping_criterion &&
    current_iteration < max_iteration_stopping_criterion)) {
    UpdateValues();
  }
}

void UnnormalizedOTSolver::UpdateValues() {
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
  if (m_x.NormInf() > divergence_tolerance || m_y.NormInf() > divergence_tolerance) {
    throw UnnormalizedOTSolverDivergenceError(*this);
  }
}

void UnnormalizedOTSolver::UpdateM() {
  double u_average;
  for (size_t k = 0; k < n_t; ++k) {
    for (size_t i = 1; i < n_x; ++i) {
      for (size_t j = 0; j < n_y; ++j) {
        u_average = (u(k, i, j) + u(k, i - 1, j)) / 2.0;
        m_x_tilde(k, i, j) = -m_x(k, i, j);
        m_x(k, i, j) = u_average / (u_average + tau_1) * (tau_1 * DxPhi(k, i, j) + m_x(k, i, j));
        m_x_tilde(k, i, j) += 2 * m_x(k, i, j);
      }
    }

    for (size_t i = 0; i < n_x; ++i) {
      for (size_t j = 1; j < n_y; ++j) {
        u_average = (u(k, i, j) + u(k, i, j - 1)) / 2.0;
        m_y_tilde(k, i, j) = -m_y(k, i, j);
        m_y(k, i, j) = u_average / (u_average + tau_1) * (tau_1 * DyPhi(k, i, j) + m_y(k, i, j));
        m_y_tilde(k, i, j) += 2 * m_y(k, i, j);
      }
    }
  }
}

void UnnormalizedOTSolver::UpdateU() {
  double a, b, c, d;
  double m_x_mean, m_y_mean;
  for (size_t k = 1; k < n_t - 1; ++k) {
    for (size_t i = 0; i < n_x; ++i) {
      for (size_t j = 0; j < n_y; ++j) {
        m_x_mean = (m_x(k, i + 1, j) + m_x(k, i, j)) / 2.0;
        m_y_mean = (m_y(k, i, j + 1) + m_y(k, i, j)) / 2.0;
       a = 1;
       b = -(tau_1 * DtPhi(k, i, j) + u(k, i, j));
       c = 0;
       d = -tau_1 / 2.0 * (pow(m_x_mean, 2) + pow(m_y_mean, 2));

        u_tilde(k, i, j) = -u(k, i, j);

        u(k, i, j) = GetRealRootsOfCubic(a, b, c, d);
        u(k, i, j) = fmax(u(k, i, j), 0);

        u_tilde(k, i, j) += 2 * u(k, i, j);
      }
    }
  }
}

void UnnormalizedOTSolver::UpdateF() {
  for (size_t k = 1; k < n_t - 1; ++k) {
    f_tilde(k) = -f(k);
    f(k) = alpha / (alpha + tau_1) *(tau_1 * GetSpacialIntegral(phi, k) + f(k));
    f_tilde(k) += 2 * f(k);
  }
}

void UnnormalizedOTSolver::UpdatePhi() {
  for (size_t k = 0; k < n_t; ++k) {
    for (size_t i = 0; i < n_x; ++i) {
      for (size_t j = 0; j < n_y; ++j) {
        phi(k, i, j) += tau_2 * (DtUTilde(k, i, j) + DivMTilde(k, i, j) - f_tilde(k));
      }
    }
  }
}

void UnnormalizedOTSolver::NormalizeU() {
  double sum_rho_0 = 0;
  for (size_t i = 0; i < n_x; ++i) {
    for (size_t j = 0; j < n_y; ++j) {
      sum_rho_0 += u(0, i, j);
    }
  }

  for (size_t k = 1; k < n_t - 1; ++k) {
    double sum_u = 0, sum_u_tilde = 0;
    for (size_t i = 0; i < n_x; ++i) {
      for (size_t j = 0; j < n_y; ++j) {
        sum_u += u(k, i, j);
        sum_u_tilde += u_tilde(k, i, j);
      }
    }
    if (sum_u > 1e-1) {
      for (size_t i = 0; i < n_x; ++i) {
        for (size_t j = 0; j < n_y; ++j) {
          u(k, i, j) *= sum_rho_0 / sum_u;
          u_tilde(k, i, j) *= sum_rho_0 / sum_u_tilde;
        }
      }
    }
  }
}

double UnnormalizedOTSolver::CheckConstraintViolation() const {
  double violation_overall = 0;
  double violation_at_time = 0;
  DoubleArray2D div_m(n_x, n_y);
  for (size_t k = 1; k < n_t - 1; ++k) {
    for (size_t i = 0; i < n_x; ++i) {
      for (size_t j = 0; j < n_y; ++j) {
        violation_at_time += fabs(DtU(k, i, j) + DivM(k, i, j) - f(k));
      }
    }
    violation_at_time *= GetXCellWidth() * GetYCellWidth() * GetTCellWidth();
    violation_overall += violation_at_time;
  }
  return violation_overall;
}

double UnnormalizedOTSolver::CheckObjectiveValue() const {
  double objective = 0;
  for (size_t k = 1; k < n_t - 1; ++k) {
    for (size_t i = 0; i < n_x; ++i) {
      for (size_t j = 0; j < n_y; ++j) {
        if (fabs(pow(m_x(k, i, j), 2) + pow(m_y(k, i, j), 2)) < 1e-8)
          continue;
        else
          objective += (pow(m_x(k, i, j), 2) + pow(m_y(k, i, j), 2)) / (u(k, i, j))
            * GetXCellWidth() * GetYCellWidth() * GetTCellWidth();
      }
    }
    objective += pow(f(k), 2) * GetTCellWidth();
  }
  return objective;
}

//double UnnormalizedOTSolver::GetDerivStencilDualityGapSpaceX() const {
//  DoubleArray3D a(1, n_x, 1), b(1, n_x, 1);
//  srand(n_x);
//  for (size_t i = 0; i < n_x; ++i) {
//    a(0, i, 0) = double(rand()) / RAND_MAX;
//    b(0, i, 0) = double(rand()) / RAND_MAX;
//  }
//
//  double integral_of_da_b = 0;
//  for (size_t i = 0; i < n_x; ++i) {
//    integral_of_da_b += SpaceDerivPrimalX(a, 0, i, 0) * b(0, i, 0);
//  }
//  integral_of_da_b *= GetXCellWidth();
//
//  double integral_of_a_db = 0;
//  for (size_t i = 0; i < n_x; ++i) {
//    integral_of_a_db -= a(0, i, 0) * SpaceDerivDualX(b, 0, i, 0);
//  }
//  integral_of_a_db *= GetXCellWidth();
//  return fabs(integral_of_a_db - integral_of_da_b);
//}
//
//double UnnormalizedOTSolver::GetDerivStencilDualityGapSpaceY() const {
//  DoubleArray3D a(1, 1, n_y), b(1, 1, n_y);
//  srand(n_y);
//  for (size_t j = 0; j < n_y; ++j) {
//    a(0, 0, j) = double(rand()) / RAND_MAX;
//    b(0, 0, j) = double(rand()) / RAND_MAX;
//  }
//
//  double integral_of_da_b = 0;
//  for (size_t j = 0; j < n_y; ++j) {
//    integral_of_da_b += SpaceDerivPrimalY(a, 0, 0, j) * b(0, 0, j);
//  }
//  integral_of_da_b *= GetYCellWidth();
//
//  double integral_of_a_db = 0;
//  for (size_t j = 0; j < n_y; ++j) {
//    integral_of_a_db -= a(0, 0, j) * SpaceDerivDualY(b, 0, 0, j);
//  }
//  integral_of_a_db *= GetYCellWidth();
//  return fabs(integral_of_a_db - integral_of_da_b);
//}

double UnnormalizedOTSolver::GetDerivStencilDualityGapTime() const {
  DoubleArray3D a(n_t, 1, 1), b(n_t, 1, 1);
  srand(n_t);
  for (size_t k = 0; k < n_t; ++k) {
    a(k, 0, 0) = double(rand()) / RAND_MAX;
    b(k, 0, 0) = double(rand()) / RAND_MAX;
  }

  double integral_of_da_b = 0;
  for (size_t k = 0; k < n_t; ++k) {
    integral_of_da_b += TimeDerivPrimal(a, k, 0, 0) * b(k, 0, 0);
  }
  integral_of_da_b *= GetTCellWidth();

  double integral_of_a_db = 0;
  for (size_t k = 0; k < n_t; ++k) {
    integral_of_a_db -= a(k, 0, 0) * TimeDerivDual(b, k, 0, 0);
  }
  integral_of_a_db *= GetTCellWidth();
  return fabs(integral_of_a_db - integral_of_da_b);
}

double UnnormalizedOTSolver::GetGradDivStencilDualityGap() {
  //DoubleArray3D a(1, n_x, n_y), b_x(1, n_x, n_y), b_y(1, n_x, n_y);
  srand(n_x*n_y);
  for (size_t i = 0; i < n_x; ++i) {
    for (size_t j = 0; j < n_y; ++j) {
      m_x(0, i, j) = double(rand()) / RAND_MAX;
      m_y(0, i, j) = double(rand()) / RAND_MAX;
      phi(0, i, j) = double(rand()) / RAND_MAX;
    }
  }

  for (size_t i = 0; i < n_x; ++i) {
    m_x(0, i, 0) = m_y(0, i, 0) = 0.0;
    m_x(0, i, n_y - 1) = m_y(0, i, n_y) = 0.0;
  }

  for (size_t j = 0; j < n_y; ++j) {
    m_x(0, 0, j) = m_y(0, 0, j) = 0.0;
    m_x(0, n_x, j) = m_y(0, n_x - 1, j) = 0.0;
  }


  double integral_of_grad_a_b = 0;
  for (size_t i = 0; i < n_x; ++i) {
    for (size_t j = 0; j < n_y; ++j) {
      integral_of_grad_a_b += phi(0, i, j) * DivM(0, i, j);
    }
  }
  integral_of_grad_a_b *= GetXCellWidth() * GetYCellWidth();
  double integral_of_a_div_b = 0;
  for (size_t i = 1; i < n_x; ++i) {
    for (size_t j = 1; j < n_y; ++j) {
      integral_of_a_div_b -= DxPhi(0, i, j) * m_x(0, i, j) + DyPhi(0, i, j) * m_y(0, i, j);
    }
  }
  for (size_t i = 0; i < n_x; ++i) {
    for (size_t j = 0; j < n_y; ++j) {
      m_x(0, i, j) = m_y(0, i, j) = phi(0, i, j) = 0.0;
    }
  }

  integral_of_a_div_b *= GetXCellWidth() * GetYCellWidth();
  return fabs(integral_of_grad_a_b - integral_of_a_div_b);
}

double UnnormalizedOTSolver::GetDivMIntegral() const {
  DoubleArray3D div_m(n_t, n_x, n_y);
  for (size_t k = 0; k < n_t; ++k) {
    for (size_t i = 0; i < n_x; ++i) {
      for (size_t j = 0; j < n_y; ++j) {
        div_m(k, i, j) = DivM(k, i, j);
      }
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
inline double UnnormalizedOTSolver::DtPhi(size_t k, size_t i, size_t j) const {
  return TimeDerivDual(phi, k, i, j);
}

inline double UnnormalizedOTSolver::DxPhi(size_t k, size_t i, size_t j) const {
  return SpaceDerivDualX(phi, k, i, j);
}

inline double UnnormalizedOTSolver::DyPhi(size_t k, size_t i, size_t j) const {
  return SpaceDerivDualY(phi, k, i, j);
}

inline double UnnormalizedOTSolver::DtU(size_t k, size_t i, size_t j) const {
  return TimeDerivPrimal(u, k, i, j);
}

inline double UnnormalizedOTSolver::DtUTilde(size_t k, size_t i, size_t j) const {
  return TimeDerivPrimal(u_tilde, k, i, j);
}

inline double UnnormalizedOTSolver::DivM(size_t k, size_t i, size_t j) const {
  return SpaceDerivPrimalX(m_x, k, i, j) + SpaceDerivPrimalY(m_y, k, i, j);
}

inline double UnnormalizedOTSolver::DivMTilde(size_t k, size_t i, size_t j) const {
  return SpaceDerivPrimalX(m_x_tilde, k, i, j) + SpaceDerivPrimalY(m_y_tilde, k, i, j);
}

inline double UnnormalizedOTSolver::GetSpacialIntegral(const DoubleArray3D & ar, size_t k) const {
  double r = 0;
  for (size_t i = 0; i < ar.GetIndex2Size(); ++i) {
    for (size_t j = 0; j < ar.GetIndex3Size(); ++j) {
      r += ar(k, i, j);
    }
  }
  return r*GetXCellWidth()*GetYCellWidth();
}

// Translations

inline double UnnormalizedOTSolver::TimeDerivPrimal(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  // return TimeDerivForwardsPrimal(ar, k, i, j);
  return TimeDerivCentralPrimal(ar, k, i, j);
  // return TimeDerivBackwardsPrimal(ar, k, i, j);
}

inline double UnnormalizedOTSolver::TimeDerivDual(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  // return TimeDerivForwardsDual(ar, k, i, j);
  return TimeDerivCentralDual(ar, k, i, j);
  // return TimeDerivBackwardsDual(ar, k, i, j);
}

inline double UnnormalizedOTSolver::SpaceDerivPrimalX(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  assert(i < n_x + 1);
  return (ar(k, i + 1, j) - ar(k, i, j)) / GetXCellWidth();
  // return SpaceDerivForwardsPrimalX(ar, k, i, j);
  // return SpaceDerivCentralPrimalX(ar, k, i, j);
  // return SpaceDerivBackwardsPrimalX(ar, k, i, j);
}

inline double UnnormalizedOTSolver::SpaceDerivDualX(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  assert(i > 0);
  return (ar(k, i, j) - ar(k, i - 1, j)) / GetXCellWidth();
  // return SpaceDerivForwardsDualX(ar, k, i, j);
  // return SpaceDerivCentralDualX(ar, k, i, j);
  // return SpaceDerivBackwardsDualX(ar, k, i, j);
}

inline double UnnormalizedOTSolver::SpaceDerivPrimalY(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  assert(j < n_y + 1);
  return (ar(k, i, j + 1) - ar(k, i, j)) / GetYCellWidth();
  // return SpaceDerivForwardsPrimalY(ar, k, i, j);
  // return SpaceDerivCentralPrimalY(ar, k, i, j);
  // return SpaceDerivBackwardsPrimalY(ar, k, i, j);
}

inline double UnnormalizedOTSolver::SpaceDerivDualY(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  assert(j > 0);
  return (ar(k, i, j) - ar(k, i, j - 1)) / GetXCellWidth();
  // return SpaceDerivForwardsDualY(ar, k, i, j);
  // return SpaceDerivCentralDualY(ar, k, i, j);
  // return SpaceDerivBackwardsDualY(ar, k, i, j);
}

// Numerical

inline double UnnormalizedOTSolver::TimeDerivCentralPrimal(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  if (k > 0 && k < n_t - 1)
    return (ar(k + 1, i, j) - ar(k - 1, i, j)) / (2 * GetTCellWidth());
  else if (k == 0)
    return (ar(1, i, j) - ar(0, i, j)) / GetTCellWidth();
  else
    return (ar(n_t - 1, i, j) - ar(n_t - 2, i, j)) / GetTCellWidth();
}

inline double UnnormalizedOTSolver::TimeDerivCentralDual(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  if (k > 1 && k < n_t - 2)
    return (ar(k + 1, i, j) - ar(k - 1, i, j)) / (2 * GetTCellWidth());
  else if (k == 0)
    return (ar(1, i, j) / 2 + ar(0, i, j)) / GetTCellWidth();
  else if (k == 1)
    return (ar(2, i, j) / 2 - ar(0, i, j)) / GetTCellWidth();
  else if (k == n_t - 2)
    return (ar(n_t - 1, i, j) - ar(n_t - 3, i, j) / 2) / GetTCellWidth();
  else
    return (-ar(n_t - 1, i, j) - ar(n_t - 2, i, j) / 2) / GetTCellWidth();
}

inline double UnnormalizedOTSolver::TimeDerivForwardsPrimal(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  if (k == n_t - 1)
    return (-ar(n_t - 1, i, j)) / GetTCellWidth();
  else
    return (ar(k + 1, i, j) - ar(k, i, j)) / GetTCellWidth();
}

inline double UnnormalizedOTSolver::TimeDerivForwardsDual(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  if (k == 0)
    return ar(0, i, j) / GetTCellWidth();
  else
    return (ar(k, i, j) - ar(k - 1, i, j)) / GetTCellWidth();
}


inline double UnnormalizedOTSolver::SpaceDerivCentralPrimalX(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  if (i > 0 && i < n_x - 1)
    return (ar(k, i + 1, j) - ar(k, i - 1, j)) / (2 * GetXCellWidth());
  if (i == 0)
    return (ar(k, 1, j) - ar(k, 0, j)) / (2 * GetXCellWidth());
  else
    return (ar(k, n_x - 1, j) - ar(k, n_x - 2, j)) / (2 * GetXCellWidth());
}

inline double UnnormalizedOTSolver::SpaceDerivCentralDualX(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  if (i > 0 && i < n_x - 1)
    return (ar(k, i + 1, j) - ar(k, i - 1, j)) / (2 * GetXCellWidth());
  if (i == 0)
    return (ar(k, 1, j) + ar(k, 0, j)) / (2 * GetXCellWidth());
  else
    return (-ar(k, n_x - 1, j) - ar(k, n_x - 2, j)) / (2 * GetXCellWidth());
}

inline double UnnormalizedOTSolver::SpaceDerivCentralPrimalY(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  if (j > 0 && j < n_y - 1)
    return (ar(k, i, j + 1) - ar(k, i, j - 1)) / (2 * GetYCellWidth());
  if (j == 0)
    return (ar(k, i, 1) - ar(k, i, 0)) / (2 * GetYCellWidth());
  else
    return (ar(k, i, n_y - 1) - ar(k, i, n_y - 2)) / (2 * GetYCellWidth());
}

inline double UnnormalizedOTSolver::SpaceDerivCentralDualY(const DoubleArray3D& ar, size_t k, size_t i, size_t j) const {
  if (j > 0 && j < n_y - 1)
    return (ar(k, i, j + 1) - ar(k, i, j - 1)) / (2 * GetYCellWidth());
  if (j == 0)
    return (ar(k, i, 1) + ar(k, i, 0)) / (2 * GetYCellWidth());
  else
    return (-ar(k, i, n_y - 1) - ar(k, i, n_y - 2)) / (2 * GetYCellWidth());
}



inline double UnnormalizedOTSolver::GetTCellWidth() const {
  return 1.0 / (n_t - 1.0);
}

inline double UnnormalizedOTSolver::GetXCellWidth() const {
  return 1.0 / (n_x - 1.0);
}

inline double UnnormalizedOTSolver::GetYCellWidth() const {
  return 1.0 / (n_y - 1.0);
}

void UnnormalizedOTSolver::SetInputs(const DoubleArray2D& mu_1, const DoubleArray2D& mu_2) {
  assert(mu_1.GetIndex1Size() == mu_2.GetIndex1Size());
  assert(mu_1.GetIndex2Size() == mu_2.GetIndex2Size());
  assert(mu_1.GetIndex1Size() == n_x);
  assert(mu_1.GetIndex2Size() == n_y);

  m_x.SetToValue(0.0);
  m_y.SetToValue(0.0);
  phi.SetToValue(0.0);
  u.SetToValue(0.0);
  f.SetToValue(0.0);

  size_t last_ind = n_t - 1;
  for (size_t i = 0; i < mu_1.GetIndex1Size(); ++i) {
    for (size_t j = 0; j < mu_1.GetIndex2Size(); ++j) {
      u(0, i, j) = u_tilde(0, i, j) = mu_1(i, j);
      u(last_ind, i, j) = u_tilde(last_ind, i, j) = mu_2(i, j);
    }
  }
}

std::ostream& operator<<(std::ostream& ostr, const UnnormalizedOTSolver& solver) {
  std::cout << "m_x = " << '\n' << solver.m_x << '\n';
  std::cout << "m_y = " << '\n' << solver.m_y << '\n';

  std::cout << "phi = " << '\n' << solver.phi << '\n';
  std::cout << "u = " << '\n' << solver.u << '\n';
  std::cout << "u_tilde = " << '\n' << solver.u_tilde << '\n';
  std::cout << "f = " << '\n' << solver.f << '\n';

  std::cout << "tau_1 = " << solver.tau_1 << '\n';
  std::cout << "tau_2 = " << solver.tau_2 << '\n';
  // std::cout << "beta = " << solver.beta << '\n';
  return ostr;
}

UnnormalizedOTSolverDivergenceError::UnnormalizedOTSolverDivergenceError(const UnnormalizedOTSolver& solver)
  : logic_error("Iterative Algorithm has diverged.") {
  iteration = solver.current_iteration;
  objective_value = solver.current_objective_value;
  constraint_value = solver.current_constraint_violation;
}

const char* UnnormalizedOTSolverDivergenceError::what() const {
  return logic_error::what();
}