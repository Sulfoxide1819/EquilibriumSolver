#include <string> 
#include <vector>
#include <Eigen/Dense>
struct Component {
  std::string name;
  double molar_mass;
  double dissociation_energy;//[J/mol]
  double vibrational_freq;
  double rotational_const;// []
  int symmetry_factor;
  bool is_atomic;
};

struct Element {
  std::string name;
  Eigen::VectorXi stoichiometry;
};

class Mixture {
public:
  Mixture(const std::vector<Component>& components, 
          const std::vector<Element> elements elements)
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
  const std::vector<Component>& components() const { return components;}
  const std::vector<Elements>& elements() const { return elements;}
private:
  std::vector<Component> components;
  std::vector<Element> elements;
  Eigen::MatrixXi stoichiometry_matrix;

};
