#include <numbers>
#include <cmath>
#include "thermo/statSum.hpp"
#include "core/constants.hpp"
#include "core/types.hpp"

#include <iostream>

double fromJANAF(const Component& comp, double T, int n);
using std::numbers::pi;
double StatSum::translational (const Component& comp, double T) {
  double base = 2 * pi * K * T / (h * h) * comp.molar_mass / NA; 
  double Z_tr = std::pow(base, 1.5);
  return std::pow(base, 1.5);
}
double StatSum::rotational(const Component& comp, double T) {
  if (comp.is_atomic) return 1.0;
  double Z_rot = K * T / (c * h * comp.symmetry_factor * comp.rotational_const);
  return Z_rot;
}

double StatSum::vibrational(const Component& comp, double T) {
  if(comp.is_atomic) return 1.0;
  double Z_vib = comp.vib_model->compute(T);
  return Z_vib;
}

double StatSum::electronic(const Component& comp, double T){
  if(comp.energy_levels.empty()){ return 1.0; }
  double Z_el = 0;
  for(const auto& level : comp.energy_levels){
    double to_exp_term = -h * c * level.energy * 100.0 / (K * T);
    if (SIMPLE_EL_STATSUM) { return level.degeneracy; }
    Z_el += level.degeneracy * std::exp(to_exp_term);
  }
  return Z_el;
}

double StatSum::total(const Component& comp, double T) {
  double exp_term = std::exp(comp.dissociation_energy / (K * T * NA));
  double Z_total = translational(comp, T) * rotational(comp, T) * vibrational(comp, T) * electronic(comp, T) * exp_term;
  return Z_total;
}

