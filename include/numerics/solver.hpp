#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include "core/types.hpp"
#include "core/constants.hpp"

namespace EquilibriumSolver {

struct SolverParameters {
  double pressure = p_0;
  double temperature;

  Eigen::VectorXd initial_mole_fractions;

  int max_iter = 0;
  double residual_tolerance = 0;
  double step_tolerance = 0;

  void validate() const {
    if (pressure <= 0) throw std::invalid_argument("Pressure must be positive");
    if (temperature <= 0) throw std::invalid_argument("Temperature must be positive");
    if (initial_mole_fractions.size() == 0)
        throw std::invalid_argument("Initial mole fractions must be provided");
    if (std::abs(initial_mole_fractions.sum() - 1.0) > 1e-10)
        throw std::invalid_argument("Initial mole fractions must sum to 1");
  }
};

struct SolverResult {
  bool success = false;
  std::string err_msg;
  int iterations = 0;
  double final_residual = 0.0;

  Eigen::VectorXd mole_fractions;
  Eigen::VectorXd concentrations;
  Eigen::VectorXd chemical_potentials;
  double log_total_density;

  double total_density() const { return /*std::exp(log_total_density);*/ concentrations.sum(); }
  double total_pressure(double T) const { return concentrations.sum() * K * T; }
};

class StatSumCache {
public:
  StatSumCache(const Mixture& mixture);

  void update(double temperature);

  const Eigen::VectorXd& get_lnZ() const { return lnZ_; }
  const Eigen::VectorXd& get_Z() const { return Z_; }
private:
  const Mixture& mixture;
  double cached_temperature = 0.0;
  Eigen::VectorXd Z_;
  Eigen::VectorXd lnZ_;
};

class EquilibriumSystem {
public: 
  EquilibriumSystem(const Mixture& mixture, 
                    const StatSumCache& statsums, 
                    const SolverParameters& params);
  
  Eigen::VectorXd compute_concentrations(const Eigen::VectorXd& x) const;
  Eigen::VectorXd compute_residuals(const Eigen::VectorXd& x, 
		                    const Eigen::VectorXd& concentrations) const;
  Eigen::MatrixXd compute_jacobian(const Eigen::VectorXd& x, 
		                   const Eigen::VectorXd& concentrations) const;

  inline int system_size() const { return Ne_ + 1; }
  const Mixture& get_mixture() const { return mixture; }
  const SolverParameters& get_params() const { return params; }
private:
  const Mixture& mixture;
  const StatSumCache& statsums;
  const SolverParameters& params;

  int Ns_;
  int Ne_;
  //Eigen::VectorXd initial_mole; //mole fractions of components [dimensionless]
  
};


class InitialGuessFinder {
  public:
    static Eigen::VectorXd find(const Mixture& mixture,
                                const StatSumCache& statsums,
                                const SolverParameters& params);    
//private:
  static Eigen::VectorXd solve_for_gamma(const Mixture& mixture,
                                          const Eigen::VectorXd& lnZ,
                                          const Eigen::VectorXd& initial_chi);
  static std::vector<int> select_equations(const Mixture& mixture,
                                          const Eigen::VectorXd& initial_chi);
};

class NewtonSolver {
public:
  struct Options {
    size_t max_iter = 1e2;
    double residual_tolerance = 1e+10;
    double step_tolerance = 1e-12;
  };

  SolverResult solve(const EquilibriumSystem& system, 
		     const Eigen::VectorXd& initial_guess);
private:
  NewtonSolver::Options options;

  bool check_convergence(const Eigen::VectorXd& residuals,
                         const Eigen::VectorXd& step,
			 const size_t& iter) const;

};

class EquilibriumCalculator {
public:
  EquilibriumCalculator(const Mixture& mixture);

  SolverResult calculate(double pressure, double temperature);
  SolverResult calculate(const SolverParameters& params);

private:
  const Mixture& mixture;
  StatSumCache statsum_cache;
};


}//namespace EquilibriumSolver


namespace EquilibriumUtils {

Eigen::VectorXd potentials_to_concentrations(const Eigen::VectorXd& gamma,
		                             const Eigen::VectorXd& lnZ,
					     const Eigen::MatrixXi& stoichiometry);

Eigen::VectorXd concentrations_to_mole_fractions(const Eigen::VectorXd& concentrations,  double total_density);

bool check_element_conservation(const Eigen::VectorXd& concentrations,
		                const Eigen::VectorXd& initial_concentrations,
                                const Eigen::MatrixXi& stoichiometry,
                                double tolerance = 1e-6);

void validate_mixture(const Mixture& mixture);
void validate_params(const EquilibriumSolver::SolverParameters& params);


}//namespace EquilibriumUtils






