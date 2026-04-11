#include <Eigen/Dense>
#include <string>
#include <cmath>
#include <algorithm>
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
    statsums(k) = std::log(StatSum::total(comp, temperature));
    k++;
  }
  this->lnZ_ = statsums;
  this->cached_temperature = temperature;
}

//EquilibriumSystem
EquilibriumSystem::EquilibriumSystem(const Mixture& mixture,
                                     const StatSumCache& statsums,
                                     const MixtureParameters& params)
  : mixture(mixture),
    statsums(statsums),
    params(params),
    Ns_(static_cast<int>(mixture.get_components().size())),
    Ne_(static_cast<int>(mixture.get_elements().size())) 
{}

Eigen::VectorXd EquilibriumSystem::compute_concentrations(const Eigen::VectorXd& gamma) const {
  double coeff = params.pressure / (K * params.temperature);
  Eigen::VectorXd to_exp =( this->statsums.get_lnZ() + 
                           this->mixture.get_stoichiometry().cast<double>()  * gamma).array();// + std::log(coeff);
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

Eigen::MatrixXd EquilibriumSystem::compute_jacobian(const Eigen::VectorXd& x,
                                                    const Eigen::VectorXd& concentrations) const {
  auto& phi = this->mixture.get_stoichiometry();
  int n = this->system_size();
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n,n);
  for (int i = 0; i < Ne_; ++i) {
    for (int j = 0; j < Ne_; ++j) {
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
    dx = jacobian.partialPivLu().solve(residuals);
    if(dx(n - 1) > 2.0) {dx(n - 1) = 2.0;}
    if(dx(n - 1) < -2.0) {dx(n - 1) = -2.0;}
    bool find_step = false;
    double lambda = 1;
    for(size_t i = 0; i < 100; ++i){
      Eigen::VectorXd x_test = x - lambda * dx;
      Eigen::VectorXd conc_test = system.compute_concentrations(x_test.head(n - 1));
      Eigen::VectorXd res_test = system.compute_residuals(x_test, conc_test);
      if(ArmijosCond(x_test, residuals, res_test, lambda)){
        find_step = true;
        x = x_test;
        break;
      } else {
        lambda *= 0.1;
      }
    }
    if(!find_step && !check_convergence(residuals, dx, iter)){throw std::runtime_error("Convergence error: unstable solution");}
    ++iter;
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

bool NewtonSolver::ArmijosCond(const Eigen::VectorXd& x_test, const Eigen::VectorXd& residuals, const Eigen::VectorXd& res_test, const double& lambda){
  static double alpha = 1e-2;
  return res_test.norm() <= (1 - lambda * alpha) * residuals.norm();
}

//InitialGuessFinder
Eigen::VectorXd InitialGuessFinder::find(const Mixture& mixture,
                                         const StatSumCache& statsums,
                                         const MixtureParameters& params) {
  int Ne = mixture.get_elements().size();
  Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(Ne + 1);
  double n_sigma_0 = params.pressure / (K * params.temperature);
  initial_guess.head(Ne) = solve_for_gamma(mixture, statsums.get_lnZ(), params.initial_mole_fractions, n_sigma_0);

  initial_guess(Ne) = std::log(n_sigma_0);
  return initial_guess;
}

Eigen::VectorXd InitialGuessFinder::solve_for_gamma(const Mixture& mixture,
                                                    const Eigen::VectorXd& lnZ,
                                                    const Eigen::VectorXd& initial_chi,
                                                    const double& n_sigma_0){
  int Ns = mixture.get_components().size();
  const Eigen::MatrixXi& phi = mixture.get_stoichiometry();
 // std::vector<int> select = select_equations(mixture, initial_chi);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(Ns);
  for(int i = 0; i < Ns; ++i){
    if(initial_chi[i] == 0) {
      b(i) = std::log(1e-9) - lnZ(i) + std::log( n_sigma_0 );
    } else {
      b(i) = std::log(initial_chi(i)) - lnZ(i) + std::log( n_sigma_0 );
    }
  }
  //Eigen::MatrixXd A = Eigen::MatrixXd::Zero(Ne, Ne);
  /*for(int i = 0; i < Ne; ++i){
    A.row(i) = phi.row(select[i]).cast<double>();
  }*/
  return phi.cast<double>().bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

}

std::vector<int> InitialGuessFinder::select_equations(const Mixture& mixture,
                                 const Eigen::VectorXd& initial_chi) {
  int n = mixture.get_components().size();
  int r = mixture.get_elements().size();
  const Eigen::MatrixXi& stoichiometry = mixture.get_stoichiometry();
  Eigen::VectorXi select = Eigen::VectorXi::Zero(n);
  select.head(r) = Eigen::VectorXi::Ones(r);

  Eigen::VectorXi indices = Eigen::VectorXi::LinSpaced(n, 0, n - 1);
  std::sort(indices.data(), indices.data() + indices.size(),
            [&initial_chi](int i1, int i2) { return initial_chi[i1] < initial_chi[i2]; });

  for(int i = r; i < n; ++i){
    double chi = initial_chi(i);
    if (chi == 0.0) continue;
    for(const auto& k : indices) {
      if (k > r - 1) continue;
      if (static_cast<bool>(select(k)) && chi > initial_chi(k) && stoichiometry(i, k) != 0) {
        select(k) = false;
        select(i) = true;
        break;
      }
    }
  }
  std::vector<int> res;
  for(int i = 0; i < select.size(); ++i){
    if(select(i) == 1) {res.push_back(i);}
  }
  return res;
}

EquilibriumCalculator::EquilibriumCalculator(const Mixture& mixture)
  : mixture(mixture), statsum_cache(mixture) {}

//EquilibriumCalculator

std::vector<SolverResult> EquilibriumCalculator::calculate(const SolverParameters& params) {
  std::vector<SolverResult> results;
  MixtureParameters mix_params;
  mix_params.pressure = params.pressure;
  mix_params.initial_mole_fractions = params.initial_mole_fractions;
  for (double T : params.temperature){
    double prev_T = this->statsum_cache.get_cached_temperature();
    mix_params.temperature = T;
    this->statsum_cache.update(T);
    EquilibriumSystem system(this->mixture, this->statsum_cache, mix_params);

    NewtonSolver solver;
    solver.set_options(params);
    
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(system.get_Ne() + 1);
    if(prev_T != 0.0 && (prev_T - T) <= 1000) {
      initial_guess.head(system.get_Ne()) = results[(results.size() - 1)].chemical_potentials;
      initial_guess(system.get_Ne()) = results[(results.size() - 1)].log_total_density;
    } else {
      initial_guess = InitialGuessFinder::find(mixture, this->statsum_cache, mix_params);
    }
    auto result = solver.solve(system, initial_guess);
  /*results.concentrations = EquilibriumUtils::potentials_to_concentrations(results.chemical_potentials,
                                                                         this->statsum_cache.get_lnZ(),
                                                                         mixture.get_stoichiometry());*/
    result.mole_fractions = EquilibriumUtils::concentrations_to_mole_fractions(result.concentrations, result.total_density());
    results.push_back(std::move(result));
  }
  return std::move(results);
}
//EquilibriumUtils

Eigen::VectorXd EquilibriumUtils::concentrations_to_mole_fractions(const Eigen::VectorXd& concentrations,  double total_density) {
  return concentrations / total_density;

}



