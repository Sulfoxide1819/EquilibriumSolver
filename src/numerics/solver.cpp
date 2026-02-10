#include <Eigen/Dense>
#include <string>
#include <cmath>
#include "core/constants.hpp"
#include "core/types.hpp"
#include "thermo/statSum.hpp"
#include "numerics/solver.hpp"

#include <iostream>

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
  double coeff = params.pressure / (K * NA * params.temperature);
  //std::cout << "size gamma: " << gamma.size() <<"\n";
  //std::cout << "size lnZ: " << statsums.get_lnZ().size() <<"\n";
  Eigen::VectorXd to_exp =( this->statsums.get_lnZ() + 
                           this->mixture.get_stoichiometry().cast<double>()  * gamma).array() + std::log(coeff);
  //std::cout <<"conc:" << "\n" << to_exp.array().exp() << "\n";
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
  //std::cout << "residual:" << "\n" << residuals << "\n";
  return residuals;
}

Eigen::MatrixXd EquilibriumSystem::compute_jacobian(const Eigen::VectorXd& x,
                                                    const Eigen::VectorXd& concentrations) const {
  auto& phi = this->mixture.get_stoichiometry();
  //std::cout << "size phi: " << phi.col(0).size() << "\n";
  int n = this->system_size();
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n,n);
  for (int i = 0; i < Ne_; ++i) {
    for (int j = 0; j < Ne_; ++j) {
      Eigen::VectorXd prod = phi.col(i).cwiseProduct(phi.col(j)).cast<double>();
      J(i,j) = prod.dot(concentrations);
    }
  }
  //std::cout << "good 1 " << "\n";
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
   || step.norm()      < options.step_tolerance 
   || iter             >= options.max_iter) { return true; }
  return false;
}

SolverResult NewtonSolver::solve(const EquilibriumSystem& system,
                                 const Eigen::VectorXd& initial_guess) {
  Eigen::VectorXd dx;
  Eigen::VectorXd residuals;
  int n = system.system_size();
  Eigen::VectorXd x = initial_guess;
  Eigen::VectorXd concentrations;
  int iter = 0;
  do {
    concentrations = system.compute_concentrations(x.head(n - 1));
    residuals = system.compute_residuals(x, concentrations);
    Eigen::MatrixXd jacobian = system.compute_jacobian(x, concentrations);
    //std::cout << "jacobian:" << jacobian << "\n";
    dx = jacobian.partialPivLu().solve(residuals);
    
    std::cout << "dx:" << dx << "\n";
    x -= dx;
    ++iter;
    std::cout << "x:"<< x << "\n";
  } while(!check_convergence(residuals, dx, iter));

  SolverResult res;
  res.iterations = iter;
  if (iter < options.max_iter) { res.success = true; }
  res.final_residual = residuals.norm();
  res.concentrations = concentrations;
  res.chemical_potentials = x.head(n - 1);
  res.log_total_density = x(n - 1);
  return res;
}

//InitialGuessFinder
Eigen::VectorXd InitialGuessFinder::find(const Mixture& mixture,
                                         const StatSumCache& statsums,
                                         const SolverParameters& params) {
  int Ne = mixture.get_elements().size();
  Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(Ne + 1);

  initial_guess(0) = std::log(params.initial_mole_fractions(2)/* * params.pressure / ( params.temperature * K)*/) - std::log(statsums.get_lnZ()(2));
  initial_guess(1) = std::log(params.initial_mole_fractions(3)/* * params.pressure / ( params.temperature * K)*/) - std::log(statsums.get_lnZ()(3));

  double n_sigma_0 = params.pressure / (K * params.temperature);
  initial_guess(Ne) = std::log(n_sigma_0);
  initial_guess(0) /= 2;
  initial_guess(1) /= 2;
  initial_guess(2) = 150.0;
  std::cout << "init_guess: "<< Ne << " " << initial_guess << "\n";
  return initial_guess;
}

EquilibriumCalculator::EquilibriumCalculator(const Mixture& mixture)
  : mixture(mixture), statsum_cache(mixture) {}

//EquilibriumCalculator

SolverResult EquilibriumCalculator::calculate(const SolverParameters& params) {
  double T = params.temperature;
  statsum_cache.update(T);
  std::cout <<  "statsums: " << statsum_cache.get_lnZ();
  EquilibriumSystem system(this->mixture, this->statsum_cache, params);
  NewtonSolver solver;
  auto results = solver.solve(system, InitialGuessFinder::find(mixture, this->statsum_cache, params));
  /*results.concentrations = EquilibriumUtils::potentials_to_concentrations(results.chemical_potentials,
                                                                         this->statsum_cache.get_lnZ(),
                                                                         mixture.get_stoichiometry());*/
  results.mole_fractions = EquilibriumUtils::concentrations_to_mole_fractions(results.concentrations, results.total_density());
  return results;
}
//EquilibriumUtils
Eigen::VectorXd EquilibriumUtils::potentials_to_concentrations(const Eigen::VectorXd& gamma,
                                                               const Eigen::VectorXd& lnZ,
                                                               const Eigen::MatrixXi& stoichiometry) {
  //Eigen::VectorXd proxy = lnZ;
  Eigen::VectorXd to_exp = lnZ + stoichiometry.cast<double>() * gamma;
  return to_exp.array().exp(); // * p/KT !!!!!!!!!!!!!!!!!!!
}

Eigen::VectorXd EquilibriumUtils::concentrations_to_mole_fractions(const Eigen::VectorXd& concentrations,
                                                                    double total_density) {
  return concentrations / total_density;

}



