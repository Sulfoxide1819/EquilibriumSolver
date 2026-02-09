#include <Eigen/Dense>
#include <string>
#include <cmath>
#include "core/constants.hpp"
#include "core/types.hpp"
#include "thermo/statSum.hpp"
#include "numerics/solver.hpp"

using namespace EquilibriumSolver;

//StatSumCache
StatSumCache::StatSumCache(const Mixture& mixture):mixture(mixture) {}

void StatSumCache::update(double temperature){
  if (temperature == cached_temperature) return;
  Eigen::VectorXd statsums(this->mixture.get_components().size());
  int k = 0;
  for (const auto& comp : this->mixture.get_components()) {
    statsums(k) = StatSum::total(comp, temperature);
    k++;
  }
  this->lnZ_ = statsums.array().log();
  this->cached_temperature = temperature;
}

//EquilibriumSystem
EquilibriumSystem::EquilibriumSystem(const Mixture& mixture,
                                     const StatSumCache& statsums,
                                     const SolverParameters& params)
  : mixture(mixture),
    statsums(statsums),
    params(params),
    Ns_(static_cast<int>(mixture.get_components().size())),
    Ne_(static_cast<int>(mixture.get_elements().size())) 
{}

Eigen::VectorXd EquilibriumSystem::compute_concentrations(const Eigen::VectorXd& gamma) const {
  //Eigen::VectorXd gamma_ = gamma.head(gamma.size() - 1);
  Eigen::VectorXd to_exp = this->statsums.get_lnZ() + 
                           this->mixture.get_stoichiometry().cast<double>()  * gamma;
  return to_exp.array().exp();
}

Eigen::VectorXd EquilibriumSystem::compute_residuals(const Eigen::VectorXd& x,
                                                     const Eigen::VectorXd& concentrations) const {
  int n = this->system_size();
  Eigen::VectorXd residuals = Eigen::VectorXd::Zero(n);
  //upgrade available:
  //init_frac * exp-> init_conc
  //conc - init_conc
  residuals.head(Ne_) = concentrations.transpose() * this->mixture.get_stoichiometry().cast<double>() - 
             (this->params.initial_mole_fractions.transpose() * this->mixture.get_stoichiometry().cast<double>()) * std::exp(x(Ne_));
  residuals(Ne_) = concentrations.sum() - params.pressure / (K * params.temperature);
  return residuals;
}

Eigen::VectorXd EquilibriumSystem::compute_jacobian(const Eigen::VectorXd& x,
                                                    const Eigen::VectorXd& concentrations) const {
  auto& phi = this->mixture.get_stoichiometry();
  int n = this->system_size();
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n,n);
  for (int i = 0; i < n-1; ++i) {
    for (int j = 0; j < n-1; ++j) {
      Eigen::VectorXd prod = phi.col(i).cwiseProduct(phi.col(j)).cast<double>();
      J(i,j) = prod.dot(concentrations);
    }
  }
  J.col(Ne_).head(Ne_) = -std::exp(x(Ne_)) * this->params.initial_mole_fractions.transpose() * phi.cast<double>();
  J.row(Ne_).head(Ne_) = concentrations.transpose() * phi.cast<double>();
  J(Ne_,Ne_) = 0.0;
  return J; 

}

//NewtonSolver
bool NewtonSolver::check_convergence(const Eigen::VectorXd& residuals,
                                     const Eigen::VectorXd& step,
                                     const size_t& iter) const {
  if (residuals.norm() < options.residual_tolerance 
   || step.norm()     < options.step_tolerance
   || iter             > options.max_iter) { return true; }
  return false;
}

SolverResult NewtonSolver::solve(const EquilibriumSystem& system,
                                 const Eigen::VectorXd& initial_guess) {
  Eigen::VectorXd dx;
  Eigen::VectorXd residuals;
  int n = system.system_size();
  Eigen::VectorXd x = initial_guess;
  int iter;
  do { 
    Eigen::VectorXd concentrations = system.compute_concentrations(x.head(n - 1));
    residuals = system.compute_residuals(x, concentrations);
    Eigen::MatrixXd jacobian = system.compute_jacobian(x, concentrations);
    dx = jacobian.partialPivLu().solve(residuals);
    x -= dx;
    ++iter;
  } while(!check_convergence(residuals, dx, iter));

  SolverResult res;
  res.iterations = iter;
  res.final_residual = residuals.norm();
  res.chemical_potentials = x.head(n - 1);
  res.log_total_density = x(n - 1);
  return res;
}

//InitialGuessFinder
Eigen::VectorXd InitialGuessFinder::find(const Mixture& mixture,
                                         const StatSumCache& statsums,
                                         const SolverParameters& params) {
  Eigen::VectorXd res;
  res << 1,1,1,1,1,1;
  return res;
}

EquilibriumCalculator::EquilibriumCalculator(const Mixture& mixture)
  : mixture(mixture), statsum_cache(mixture) {}

//EquilibriumCalculator

SolverResult EquilibriumCalculator::calculate(const SolverParameters& params) {
  EquilibriumSystem system(this->mixture, this->statsum_cache, params);
  NewtonSolver solver;
  auto results = solver.solve(system, InitialGuessFinder::find(mixture, this->statsum_cache, params));
  results.concentrations = EquilibriumUtils::potentials_to_concentrations(results.chemical_potentials,
                                                                         this->statsum_cache.get_lnZ(),
                                                                         mixture.get_stoichiometry());
  if (EquilibriumUtils::check_element_conservation(results.concentrations, system.get_params().initial_mole_fractions * results.total_density(), mixture.get_stoichiometry())) {
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
                                                               const Eigen::MatrixXi& stoichiometry) {
  //Eigen::VectorXd proxy = lnZ;
  Eigen::VectorXd to_exp = lnZ + stoichiometry.cast<double>() * gamma;
  return to_exp.array().exp();
}

Eigen::VectorXd EquilibriumUtils::concentrations_to_mole_fractions(const Eigen::VectorXd& concentrations,
                                                                    double total_density) {
  return concentrations / total_density;

}

bool EquilibriumUtils::check_element_conservation(const Eigen::VectorXd& concentrations,
                                                  const Eigen::VectorXd& initial_concentrations,
                                                  const Eigen::MatrixXi& stoichiometry,
                                                  double tolerance) {
  Eigen::VectorXd error = (concentrations - initial_concentrations).transpose() * stoichiometry.cast<double>();
  return error.norm() < tolerance;
}












