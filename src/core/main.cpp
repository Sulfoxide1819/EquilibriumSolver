#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <Eigen/Dense>
#include "core/constants.hpp"
#include "core/types.hpp"
#include "thermo/statSum.hpp"
#include "numerics/solver.hpp"
#include "core/parse_utils.hpp"
#include "core/input.hpp"

using namespace EquilibriumSolver;
using namespace std;

int main(int argc, char* argv[]) {
  try {
    if(argc == 0) throw std::runtime_error("First argument must be configuration file");
    ifstream config_file(argv[1]);
    nlohmann::json config = nlohmann::json::parse(config_file);
    Read read(config);
    vector<Component> components = read.components();
    cout << "=== Starting Equilibrium Solver ===" << endl;
    SolverParameters params = read.params();

    ifstream f("../data/molecules.json");
    nlohmann::json data = nlohmann::json::parse(f);
    nlohmann::json data_thermo = nlohmann::json::parse(ifstream("../data/gibbs_entalpy.json"));
    for(auto& comp : components){
        MixtureBuild::pull_properties(data, comp);
    }
    vector<Element> elements = MixtureBuild::get_elements(data, components);
    for(auto& el : elements) {
      MixtureBuild::get_stoichiometry(data, el, components);
    }   
 // Mixture init
    Mixture mixture(components, elements);
    //Calculator init
    EquilibriumCalculator calculator(mixture);

    
    
    params.max_iter = 1e4;
    params.residual_tolerance = 1e+11;
    params.step_tolerance = 1e-16;

//OUTPUT
    vector<SolverResult> results = calculator.calculate(params);
    SolverResult result = results[0];

    filesystem::path path(argv[1]);
    ofstream outfile(path.replace_extension(".csv"));
    if(!outfile.is_open()){
      throw runtime_error("Cannot open file for writing!");
    }

    outfile << "T";
    size_t size = components.size();
    for(size_t i = 0; i < size; ++i){
      outfile << "," << components[i].name;
    }
    outfile << endl;

    size_t k = 0;
    for(SolverResult res : results){
      outfile << params.temperature[k];
      for(size_t i = 0; i < size; ++i){
        outfile << "," << res.mole_fractions(i);
      }
      outfile << endl;
      ++k;
    }
    outfile.close();

    cout << "=========================================\n";
    cout << "Equilibrium Composition Calculation\n";
    cout << "=========================================\n";
    cout << "Parameters:\n";
    cout << "  Pressure: " << params.pressure << " Pa\n";
    cout << "  Temperature: " << params.temperature[0] << " K\n";
    cout << "  Initial composition (mole fractions):\n";
    for (size_t i = 0; i < components.size(); ++i) {
      cout << "    " << components[i].name << ": " << params.initial_mole_fractions(i) << "\n";
    }
    cout << "\nResults:\n";
    cout << "  Success: " << (result.success ? "Yes" : "No") << "\n";
    if (!result.success) {
      cout << "  Error: " << result.err_msg << "\n";
    }
    cout << "  Iterations: " << result.iterations << "\n";
    cout << "  Final residual: " << result.final_residual << "\n";
    cout << "  Total density: " << result.total_density() << " particles/m^3\n";
    cout << "\nEquilibrium composition:\n";
       
    for (int i = 0; i < result.mole_fractions.size(); ++i) {
      cout << "  " << components[i].name << ":\n";
      cout << "    Mole fraction: " << result.mole_fractions(i) << "\n";
      cout << "    Concentration: " << result.concentrations(i) << " particles/m^3\n";
    }
      
    cout << "\nChemical potentials (gamma):\n";
    for (int i = 0; i < result.chemical_potentials.size(); ++i) {
      cout << "  " << elements[i].name << ": " << result.chemical_potentials(i) << "\n";
    }
        
  } catch (const exception& e) {
    cerr << "Error: " << e.what() << endl;
    return 1;
  }
 
  return 0;
}
