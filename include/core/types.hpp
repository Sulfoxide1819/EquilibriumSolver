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
  Eigen::VectorXd stoichiometry;
};

class Mixture {
public:
  const std::vector<Component>& get_components { return components;}
private:
  std::vector<Component> components;
  std::vector<Element> elements;
  Eigen::MatrixXi stoichiometry_matrix;

};
