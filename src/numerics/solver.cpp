#include "constant.hpp"
#include "types.hpp"
#include "statSum.hpp"
#include "solver.hpp"

using namespace EquilibriumSolver;

//StatSumCache
StatSumCache::StatSumCache(const Mixture& mixture):mixture(mixture) {}

void StatSumCache::update(double temperature){
  Eigen::VectorXd statsums(mixture.get_size_comp());
  int k = 0;
  for (const auto comp& : this->mixture.get_components()) {
    statsums(k) = StatSum::total(comp, temperature);
    k++;
  }
  this->lnZ_ = statsums.log();
  this->cached_temperature = temperature;
}

//EquilibriumSystem
EquilibriumSystem::EquilibriumSystem(

Eigen::VectorXd EquilibriumSystem::compute_concentrations(const Eigen::VectorXd& gamma) const {
  Eigen::VectorXd gamma_ = gamma.head(gamma.size() - 1);
  Eigen::VectorXd to_exp = this->statsums.get_lnZ() + 
                           this->mixture.stoichiometry_matrix * gamma_;
  return to_exp.exp();
}

Eigen::VectorXd EquilibriumSystem::compute_residuals(const Eigen::VectorXd& gamma, 
                                                     const Eigen::VectorXd& concentrations) const {
  Eigen::VectorXd residuals = concentrations.transpose() * this->mixture.stoichiometry_matrix - 
                              (this->initial_distribution.transpose() * this->mixture.stoichiometry_matrix) * (std::exp(gamma(gamma.size()-1)));
  residuals.insert(last) = concentrations.sum() - params.pressure / (K * params.temperature);
  return residuals;
}

Eigen::VectorXd EquilibriumSystem::compute_jacobian(const Eigen::VectorXd& gamma,
                                                    const Eigen::VectorXd& concentrations) const {
  auto& m = this->mixture.stoichiometry_matrix;
  int n = this->system_size();
  Eigen::MatrixXd A(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      A(i,j) = m.col(i).cwiseProduct(m.col(j)).transpose() * concentrations;
    }
  }
  A.col(n).head(n-1) = -std::exp(gamma(n)) * this->initial_distribution.transpose() * m;
  A.row(n).head(n-1) = concentrations.transpose() * m;
  A(n,n) = 0.0;
  return A; 

}

//NewtonSolver
bool NewtonSolver::check_convergence(const Eigen::VectorXd& residuals,
                                     const Eigen::VectorXd& step,
                                     const size_t& iter) const {
  if (residuals.norm() < options.residual_tolerance 
   || step.norm()     < options.step_tolerance
   || iter             > max_iter) { return true; }
  return false;
}

SolverResult NewtonSolver::solve(const EquilibriumSystem& system,
                                 const Eigen::VectorXd& initial_guess) {
  Eigen::VectorXd step;
  Eigen::VectorXd residuals;
  Eigen::VectorXd gamma = initial_guess;
  int iter;
  do { 
    Eigen::VectorXd concentrations = system.compute_concentrations();
    residuals = system.compute_residuals(gamma, concentrations);
    Eigen::MatrixXd jacobian = system.compute_jacobian(gamma, concentration);
    step = jacobian.partialPivLu().solve(residuals);
    gamma -= step;
    ++iter;
  } while(!check_convergence(residuals, step, iter))

  SolverResult res;
  res.iterations = iter;
  res.final_residual = residuals;
  res.chemical_potentials = gamma.head(gamma.size()-2);
  res.total_density = gamma.tail(1);
  return res;
}

//EquilibriumCalculator

SolverResult EquilibriumCalculator::calculate(double pressure, double temperature) {
  EquilibriumSystem system(this->mixture, this->statsum_cache, SolverParameters(pressure, temperature));
  NewtonSolver solver();
  auto results = solver.solve(system, InitialGuessFinder.find());
  results.concentration = EquilibriumUtils::potentials_to_concentrations(results.chemical_potentials,
                                                                         statsum_cache,
                                                                         mixture.stoichiometry_matrix);
  if (EquilibriumUtils::check_element_conservation(results.concentration, system.initial_ditribution * results.total_density, mixture.stoichiometry_matrix)) {
  results.success = true;
  } else {
    results.err_msg = "No conservation of elements";
    return results;
  }
  results.mole_fractions = EquilibriumUtils::concentrations_to_mole_fractions(results.concentration, results.total_density);
  return results;
}
//EquilibriumUtils
Eigen::VectorXd EquilibriumUtils::potentials_to_concentrations(const Eigen::VectorXd& gamma,
                                                               const Eigen::VectorXd& lnZ,
                                                               const Eigen::MatrixXi& stoichiometry) const {
  Eigen::VectorXd to_exp = lnZ + stoichiometry * gamma;
  return to_exp.exp();
}

Eigen::VectorXd& EquilibriumUtils::concentrations_to_mole_fractions(const Eigen::VectorXd& concentrations,
                                                                    double total_density) const {
  return concentrations / total_density;

}

bool EquilibriumUtils::check_element_conservation(const Eigen::VectorXd& concentrations,
                                                  const Eigen::VectorXd& initial_concentrations,
                                                  const Eigen::MatrixXi& stoichiometry,
                                                  double tolerance = 1e-10) {
  auto error = (concentrations - initial_concentrations).transpose() * stoichiometry;
  return error.norm() < tolerance;
}


