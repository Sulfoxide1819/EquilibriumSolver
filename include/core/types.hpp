#pragma once
#include <string> 
#include <vector>
#include <Eigen/Dense>
struct Component {
  Component(std::string name):name(name){}
  
  std::string name;
  double molar_mass = 0;
  double dissociation_energy = 0;//[J/mol]
  std::vector<double> vibrational_freq;
  double rotational_const = 0;// []
  int symmetry_factor = 1;
  
  bool is_atomic;
  bool is_charged;

  struct EnergyLevel {
    std::string state;
    double energy;
    int degeneracy;
  };
  std::vector<EnergyLevel> energy_levels;
};


struct Element {
  Element(std::string name): name(name){}
  std::string name;
  std::vector<int> stoichiometry;
  //Eigen::VectorXi stoichiometry;
};

class Mixture {
public:
  Mixture(const std::vector<Component>& components, 
          const std::vector<Element>& elements)
  : components(components), elements(elements) {

  int n = components.size();
  int m = elements.size();
  Eigen::MatrixXi phi = Eigen::MatrixXi::Zero(n,m);
  for(int j = 0; j < m; ++j) {
    for(int i = 0; i < n; ++i) {
      phi(i,j) = elements[j].stoichiometry[i];
    }
  }
  this->stoichiometry_matrix = phi;
}
  const std::vector<Component>& get_components() const { return components;}
  const std::vector<Element>& get_elements() const { return elements;}
  const Eigen::MatrixXi& get_stoichiometry() const { return stoichiometry_matrix; }
private:
  std::vector<Component> components;
  std::vector<Element> elements;
  Eigen::MatrixXi stoichiometry_matrix;

};
