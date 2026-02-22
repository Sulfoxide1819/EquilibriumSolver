#pragma once
#include "core/types.hpp"
#include <vector>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include "numerics/solver.hpp"


class Read {
public:
  explicit Read(const nlohmann::json& config): config(config){}
  std::vector<Component> components();
  EquilibriumSolver::SolverParameters params();
private:
  const nlohmann::json& config;
};
