#pragma once
#include <Eigen/Dense>

namespace EquilibriumSolver {

struct SolverParameters {
  double pressure;
  double temperature;

  int max_iter = 0;
  double absolute_tolerance = 0.0;
  double relative_tolerance = 0.0;
};

struct SolverResults {
  bool success = false;
  std::string err_msg;
  int iterations = 0;
  double final_residual = 0.0;

  Eigen::VectorXd mole_fractions;
  Eigen::VectorXd concentrations;
  Eigen::VectorXd chemical_potentials;

  double total_density;
  double total_pressure;
};

class StatSumCache {
public:
  StatSumCache(const Mixture& mixture);

  void update(double temperature);

  const Eigen::VectorXd& get_lnZ() const { return lnZ_; }
private:
  const Mixture& mixture;
  double cached_temperature = 0.0;
  Eigen::VectorXd lnZ_;
};

class EquilibriumSystem {
public: 
  EquilibriumSystem(const Mixture& mixture, const StatSumCache& statsums, const SolverParameters& params);
  
  Eigen::VectorXd compute_concentrations(const Eigen::VectorXd& gamma) const;
  Eigen::VectorXd compute_residuals(const Eigen::VectorXd& gamma, 
		                    const Eigen::VectorXd& concentrations);
  Eigen::VectorXd compute_jacobian(const Eigen::VectorXd& gamma, 
		                   const Eigen::VectorXd& concentrations) const;
private:
  const Mixture& mixture;
  const StatSumCache& statsums;
  const SolverParameters& params;
  Eigen::VectorXd initial_distribution; //mole fractions of components [dimensionless]
};

class InitialGuessFinder {
  

};

class NewtonSolver {
public:
  struct Options {
    size_t max_iter = 100;
    double absolute_tolerance = 1e-10;
    double relative_tolerance = 1e-8;
  };

  SolverResult solve(const EquilibriumSystem& system, 
		     const Eigen::VectorXd& initial_guess);
private:
  Options options;

  bool check_convergence(const Eigen::VectorXd& residuals,
		         const Eigen::VectorXd& error,
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
					     const Eigen::MatrixXi& stoichiometry) const;

Eigen::MatrixXi& concentrations_to_mole_fractions(const Eigen::VectorXd& concentrations,
                                                  double total_density) const;

bool check_element_conservation(const Eigen::VectorXd& concentrations,
		                const Eigen::VectorXd& initial_concentrations,
                                const Eigen::MatrixXi& stoichiometry,
                                double tolerance = 1e-10);

void validate_mixture(const Mixture& mixture);
void validate_params(const SolverParameters& params);


}//namespace EquilibriumUtils






