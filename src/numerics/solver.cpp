#include <Eigen/Dense>
#include <string>
#include <cmath>
#include "constant.hpp"
#include "types.hpp"
#include "statSum.hpp"
#include "solver.hpp"

using namespace EquilibriumSolver;

//StatSumCache
StatSumCache::StatSumCache(const Mixture& mixture):mixture(mixture) {}

void StatSumCache::update(double temperature){
  if (temperature == cached_temperature) return;
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
EquilibriumSystem::EquilibriumSystem(const Mixture& mixture,
                                     const StatSumCache& statsums,
                                     const SolverParameters& params)
  : mixture(mixture),
    statsums(statsums),
    params(params),
    Ns_(static_cast<int>(mixture.components.size())),
    Ne_(static_cast<int>(mixture.elements.size())) 
{}

Eigen::VectorXd EquilibriumSystem::compute_concentrations(const Eigen::VectorXd& gamma) const {
  //Eigen::VectorXd gamma_ = gamma.head(gamma.size() - 1);
  Eigen::VectorXd to_exp = this->statsums.get_lnZ() + 
                           this->mixture.stoichiometry_matrix * gamma;
  return to_exp.exp();
}

Eigen::VectorXd EquilibriumSystem::compute_residuals(const Eigen::VectorXd& x,
                                                     const Eigen::VectorXd& concentrations) const {
  int n = this->system_size();
  Eigen::VectorXd residuals = Eigen::VectorXd::Zero(n);
  residuals.head(n) = concentrations.transpose() * this->mixture.stoichiometry_matrix - 
             (this->initial_mole_fractions.transpose() * this->mixture.stoichiometry_matrix) * std::exp(x(n));
  residuals(n) = concentrations.sum() - params.pressure / (K * params.temperature);
  return residuals;
}

Eigen::VectorXd EquilibriumSystem::compute_jacobian(const Eigen::VectorXd& x,
                                                    const Eigen::VectorXd& concentrations) const {
  auto& phi = this->mixture.stoichiometry_matrix;
  int n = this->system_size();
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n,n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      J(i,j) = phi.col(i).cwiseProduct(phi.col(j)).transpose() * concentrations;
    }
  }
  J.col(n).head(n-1) = -std::exp(x(n)) * this->initial_mole_fractions.transpose() * phi;
  J.row(n).head(n-1) = concentrations.transpose() * phi;
  J(n,n) = 0.0;
  return J; 

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
  Eigen::VectorXd dx;
  Eigen::VectorXd residuals;
  
  Eigen::VectorXd x = initial_guess;
  int iter;
  do { 
    Eigen::VectorXd concentrations = system.compute_concentrations(x.head(Ne_));
    residuals = system.compute_residuals(x, concentrations);
    Eigen::MatrixXd jacobian = system.compute_jacobian(x, concentration);
    dx = jacobian.partialPivLu().solve(residuals);
    x -= dx;
    ++iter;
  } while(!check_convergence(residuals, dx, iter))

  SolverResult res;
  res.iterations = iter;
  res.final_residual = residuals;
  res.chemical_potentials = x.head(Ne_);
  res.log_total_density = x.tail(1);
  return res;
}

//EquilibriumCalculator

SolverResult EquilibriumCalculator::calculate(const SolverParameters& params) {
  EquilibriumSystem system(this->mixture, this->statsum_cache, params);
  NewtonSolver solver;
  auto results = solver.solve(system, InitialGuessFinder.find());
  results.concentration = EquilibriumUtils::potentials_to_concentrations(results.chemical_potentials,
                                                                         statsum_cache,
                                                                         mixture.stoichiometry_matrix);
  if (EquilibriumUtils::check_element_conservation(results.concentrations, system.initial_mole_fractions * results.total_density(), mixture.stoichiometry_matrix)) {
  results.success = true;
  } else {
    results.err_msg = "No conservation of elements";
    return results;
  }
  results.mole_fractions = EquilibriumUtils::concentrations_to_mole_fractions(results.concentrations, results.total_density());
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












