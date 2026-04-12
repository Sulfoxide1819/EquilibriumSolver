#include "core/input.hpp"
#include <string>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
using EquilibriumSolver::SolverParameters;

  std::vector<Component> Read::components() {
    if (!config.contains("components")) throw std::runtime_error("No 'components' section in config");
    if (!config["components"].is_array() || config["components"].empty()) throw std::runtime_error("'components' must be a non-empty array");
    std::vector<Component> comps;
    for(std::string& name : config.at("components").get<std::vector<std::string>>()){
      comps.push_back(std::move(Component(name)));
    }
    return comps;
}

  SolverParameters Read::params() {
    /*if(!config.contains("Parameters")) throw std::runtime_error("No 'Parameters' section in config");
    if(config["parameters"].empty()) throw std::runtime_error("'parameters' must be a non-empty");*/
    SolverParameters p;
    if (!config.contains("temperature")) throw std::runtime_error("No 'temperature' in 'Parameters'");
    if (config["temperature"].empty()) throw std::runtime_error("'temperature' must be given");
    if (config["temperature"].is_array()){
      p.temperature = config.at("temperature").get<std::vector<double>>();
    } else if (config["temperature"].is_number()){
      p.temperature.push_back(config.at("temperature").get<double>());
    } else if (config["temperature"].is_string()) {
      try {
        p.temperature = parseTemperatureRange(config["temperature"].get<std::string>());
      } catch (const std::exception& e) {
        std::cerr << "Configuration error: temperature: ";
        std::cerr << e.what() << std::endl;
        exit(1);
      }
    } else throw std::runtime_error("Configuration error: temperature: format error!");
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
  std::vector<double> parseTemperatureRange(const std::string& spec) {
    std::vector<std::string> parts;
    std::stringstream ss(spec);
    std::string item;
    while (std::getline(ss, item, ':')) {
      parts.push_back(item);
    }

    if (parts.size() != 3){ throw std::runtime_error("Incorrect format of the field 'temperature'"); }

    double start = std::stod(parts[0]);
    double end = std::stod(parts[1]);
    double step = std::stod(parts[2]);

    if (start <= 0) throw std::runtime_error("Start cannot be zero");
    if (step <= 0) throw std::runtime_error("Step must be positive");
    if (start > end) throw std::runtime_error("Start must be less than end");

    std::vector<double> temps;
    for (size_t i = start; i <= end + 1e-9 * step; i += step) {
      temps.push_back(i);
    }
    return temps;
  }
