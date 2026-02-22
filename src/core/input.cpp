#include "core/input.hpp"
#include <string>
using EquilibriumSolver::SolverParameters;

  std::vector<Component> Read::components() {
    if(!config.contains("components")) throw std::runtime_error("No 'components' section in config");
    if(!config["components"].is_array() || config["components"].empty()) throw std::runtime_error("'components' must be a non-empty array");
    std::vector<Component> comps;
    for(std::string& name : config.at("components").get<std::vector<std::string>>()){
      comps.push_back(Component(name));
    }
    return comps;
}

  SolverParameters Read::params() {
    /*if(!config.contains("Parameters")) throw std::runtime_error("No 'Parameters' section in config");
    if(config["parameters"].empty()) throw std::runtime_error("'parameters' must be a non-empty");*/
    SolverParameters p;
    if(!config.contains("temperature")) throw std::runtime_error("No 'temperature' in 'Parameters'");
    if(config["temperature"].empty()) throw std::runtime_error("'temperature' must be given");
    p.temperature = config.at("temperature").get<double>();
    /* 
    if(!config.contains("pressure")) throw std::runtime_error("No 'pressure' in 'Parameters'");
    if(config["pressure"].empty()) throw std::runtime_error("'pressure' must be given");*/
    if(config.contains("pressure") && !config["pressure"].empty()) {
      p.pressure = config.at("pressure").get<double>();
    }

    if(!config.contains("initial_mole_fractions")) throw std::runtime_error("No 'initial_mole_fractions' in 'Parameters'");
    if(config["initial_mole_fractions"].empty()) throw std::runtime_error("'initial_mole_fractions' must be given");
    
    if(config["initial_mole_fractions"].size() != config["components"].size()) throw std::runtime_error("Length of 'initial_mole_fractions' must be equal to numbers of components");
    std::vector chi_v = config.at("initial_mole_fractions").get<std::vector<double>>();
    
    p.initial_mole_fractions = Eigen::Map<Eigen::VectorXd>(
      chi_v.data(), config["initial_mole_fractions"].size());
    
    return p;
   }
    
