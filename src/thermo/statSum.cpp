#include <numbers>
#include <cmath>
#include "thermo/statSum.hpp"
#include "core/constants.hpp"
#include "core/types.hpp"

using std::numbers::pi;
double StatSum::translational (const Component& comp, double T) {
  double base = 2 * pi * K * T / (h * h) * comp.molar_mass / NA; 
  double Z_tr = std::pow(base, 1.5);
  //std::cout << "Z_tr ( " << comp.name << " ): " << Z_tr << "\n";
  return std::pow(base, 1.5);
}
double StatSum::rotational(const Component& comp, double T) {
  if (comp.is_atomic) return 1.0;
  return K * T / (c * h * comp.symmetry_factor * comp.rotational_const);
}
double StatSum::vibrational(const Component& comp, double T) {
  if (comp.is_atomic) return 1.0;
  double denom = 1 - std::exp(-h * c * comp.vibrational_freq / (K * T));
  return (denom > 0) ? 1/denom : 1.0;
}

double StatSum::electronic(const Component& comp, double T){
  if(!comp.is_atomic){ return 1.0; }
  double Z_el = 0;
  for(const auto& level : comp.energy_levels){
    double to_exp_term = -h * c * level.energy / (K * T);
    Z_el += level.degeneracy * std::exp(to_exp_term);
  }
  return Z_el;
}


double StatSum::total(const Component& comp, double T) {
  double exp_term = std::exp(comp.dissociation_energy / (K * T * NA));
  return translational(comp, T) * rotational(comp, T) * vibrational(comp, T) * electronic(comp, T) * exp_term;
}
