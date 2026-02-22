#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
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
    if(argc == 0) throw std::runtime_error("First argument must be configuration file");
cout << 53348.4 * 100 * NA * h * c;
    ifstream config_file(argv[1]);
    nlohmann::json config = nlohmann::json::parse(config_file);
    Read read(config);
    vector<Component> components = read.components();
    cout << "=== Starting Equilibrium Solver ===" << endl;
    SolverParameters params = read.params();
/*
    SolverParameters params;
    params.pressure = 1 * 101.325; // Па
    params.temperature = 3000.0; // K
*/ 
    //vector<Component> components = {Component("N"),Component("O"), Component("N2"), Component("O2"), Component("NO")};
    //vector<Component> components = {Component("e-"),Component("N"),Component("O"), Component("Ar"), Component("N2"), Component("O2"),Component("NO"),Component("N+"),Component("O+"),Component("Ar+"),Component("NO+")};
    //vector<Component> components = {Component("C"),Component("O"), Component("CO"), Component("CO2"), Component("O2")};

    ifstream f("../data/molecules.json");
    nlohmann::json data = nlohmann::json::parse(f);
    nlohmann::json data_thermo = nlohmann::json::parse(ifstream("../data/gibbs_entalpy.json"));
    for(auto& comp : components){
//      if(MixtureBuild::pull_thermo(data_thermo, comp, params.temperature)){
        MixtureBuild::pull_properties(data, comp);
//      }
    }
    //vector<Element> elements = {Element("N"), Element("O")};
    //vector<Element> elements = {Element("e-"), Element("N"), Element("O"), Element("Ar")};
    //vector<Element> elements = {Element("C"), Element("O")};
    vector<Element> elements = MixtureBuild::get_elements(data, components);
    for(auto& el : elements) {
      MixtureBuild::get_stoichiometry(data, el, components);
    }   
 // Mixture init
    Mixture mixture(components, elements);
    cout << mixture.get_stoichiometry() << "\n";
    //Calculator init
    EquilibriumCalculator calculator(mixture);

    
    //Eigen::VectorXd initial_mole_frac(components.size());
    //initial_mole_frac << 0.0, 0.0, 0.8, 0.2, 0.0;
    //initial_mole_frac << 0.0, 0.0, 0.0, 0.01, 0.78, 0.21 , 0.0, 0.0, 0.0, 0.0, 0.0;
    //initial_mole_frac << 0.03, 0.45, 0.51, 0.005, 0.005; 
    //initial_mole_frac << 0.00, 0.33, 0.33, 0.33, 0.00;
    //params.initial_mole_fractions = initial_mole_frac;
    
    params.max_iter = 100;
    params.residual_tolerance = 1e-10;
    params.step_tolerance = 1e-8;
    
    try {
        SolverResult result = calculator.calculate(params);
        
        cout << "=========================================\n";
        cout << "Equilibrium Composition Calculation\n";
        cout << "=========================================\n";
        cout << "Parameters:\n";
        cout << "  Pressure: " << params.pressure << " Pa\n";
        cout << "  Temperature: " << params.temperature << " K\n";
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
