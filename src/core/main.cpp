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

using namespace EquilibriumSolver;
using namespace std;

int main() {
    cout << "=== Starting Equilibrium Solver ===" << endl;

    //vector<Component> components = {Component("N"),Component("O"), Component("N2"), Component("O2"), Component("NO")};
    //vector<Component> components = {Component("e-"),Component("N"),Component("O"), Component("Ar"), Component("N2"), Component("O2"),Component("NO"),Component("N+"),Component("O+"),Component("Ar+"),Component("NO+")};
    vector<Component> components = {Component("C"),Component("O"), Component("CO"), Component("CO2"), Component("O2")};
    // Создаем компоненты для 5-компонентного воздуха (N, O, N2, O2, NO)
    /*vector<Component> components = {
        // N (атомный азот)
        {
            "N",                    // name
            14.0067e-3,            // molar_mass [kg/mol]
            0.0,                    // dissociation_energy [J/mol]
            0.0,                    // vibrational_freq [1/m]
            0.0,                    // rotational_const [1/m]
            1,                      // symmetry_factor
            true                    // is_atomic
        },
        // O (атомный кислород)
        {
            "O",                    // name
            15.9994e-3,            // molar_mass [kg/mol]
            0.0,                    // dissociation_energy [J/mol]
            0.0,                    // vibrational_freq [1/m]
            0.0,                    // rotational_const [1/m]
            1,                      // symmetry_factor
            true                    // is_atomic
        },
        // N2
        {
            "N2",                   // name
            28.0134e-3,            // molar_mass [kg/mol]
            9.4158e5,              // dissociation_energy [J/mol] из таблицы
            2.35857e5,             // vibrational_freq [1/m]
            199.8,                 // rotational_const [1/m]
            2,                      // symmetry_factor (гомоядерный)
            false                   // is_atomic
        },
        // O2
        {
            "O2",                   // name
            31.9988e-3,            // molar_mass [kg/mol]
            4.9358e5,              // dissociation_energy [J/mol]
            1.58019e5,             // vibrational_freq [1/m]
            143.77,                // rotational_const [1/m]
            2,                      // symmetry_factor (гомоядерный)
            false                   // is_atomic
        },
        // NO
        {
            "NO",                   // name
            30.0061e-3,            // molar_mass [kg/mol]
            6.2685e5,              // dissociation_energy [J/mol]
            2.16981e5,             // vibrational_freq [1/m]
            167.20,                // rotational_const [1/m]
            1,                      // symmetry_factor (гетероядерный)
            false                   // is_atomic
        }
    };*/

    // Создаем элементы (N и O)
    // Стехиометрия: для каждого элемента вектор из 5 значений (N, O, N2, O2, NO)
    vector<Element> elements1 = {
    
            Element("N"),
            //(Eigen::VectorXi(5) << 1, 0, 2, 0, 1).finished() // N в N, O, N2, O2, NO
            //{1,0,2,0,1}
        
        
            Element("O")
            //(Eigen::VectorXi(5) << 0, 1, 0, 2, 1).finished() // O в N, O, N2, O2, NO
            //{0,1,0,2,1}
        
    };
    ifstream f("../data/molecules.json");
    nlohmann::json data = nlohmann::json::parse(f);
    for(auto& comp : components){
      MixtureBuild::pull_properties(data, comp);
    }
    //vector<Element> elements = {Element("e-"), Element("N"), Element("O"), Element("Ar")};
    vector<Element> elements = {Element("C"), Element("O")};
    for(auto& el : elements) {
      MixtureBuild::get_stoichiometry(data, el, components);
    }   
 // Создаем смесь
    Mixture mixture(components, elements);
    cout << mixture.get_stoichiometry() << "\n";
    //cout << InitialGuessFinder::select_equations(mixture, (Eigen::VectorXd(5) << 0.0, 0.0, 0.2, 0.8, 0.0).finished());
    // Создаем калькулятор
    EquilibriumCalculator calculator(mixture);

    // Параметры расчета: p = 0.001 * p0 = 101.325 Па, T = 1500 K
    SolverParameters params;
    params.pressure = 1 * 101.325; // Па
    params.temperature = 500.0; // K
    
    // Начальные мольные доли (недиссоциированный воздух: 80% N2, 20% O2)
    Eigen::VectorXd initial_mole_frac(components.size());
    //initial_mole_frac << 0.0, 0.0, 0.8, 0.2, 0.0;
    //initial_mole_frac << 0.0, 0.0, 0.0, 0.01, 0.78, 0.21 , 0.0, 0.0, 0.0, 0.0, 0.0;
    initial_mole_frac << 0.4, 0.4, 0.2, 0.0, 0.0; 
    params.initial_mole_fractions = initial_mole_frac;
    
    params.max_iter = 100;
    params.residual_tolerance = 1e-10;
    params.step_tolerance = 1e-8;
    
    try {
        // Выполняем расчет
        SolverResult result = calculator.calculate(params);
        
        // Выводим результаты
        cout << "=========================================\n";
        cout << "Equilibrium Composition Calculation\n";
        cout << "=========================================\n";
        cout << "Parameters:\n";
        cout << "  Pressure: " << params.pressure << " Pa\n";
        cout << "  Temperature: " << params.temperature << " K\n";
        cout << "  Initial composition (mole fractions):\n";
        for (size_t i = 0; i < components.size(); ++i) {
            cout << "    " << components[i].name << ": " << initial_mole_frac(i) << "\n";
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
        
        // Выводим мольные доли и концентрации
        //cout << fixed << setprecision(6);
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
