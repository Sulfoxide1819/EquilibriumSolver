#include <numbers>
#include <cmath>
#include "thermo/statSum.hpp"
#include "core/constants.hpp"
#include "core/types.hpp"

#include <iostream>

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
  const std::vector<double>& vibf = comp.vibrational_freq;
  if(vibf.size() > 1) {
    std::vector<double> vibl; //energy of vibrational levels
    double de = comp.dissociation_energy / NA;
    std::vector<int> lim = {34, 67, 20};
    for(int i = 0; i < lim[0]; ++i){
      for(int j = 0; j < lim[1]; ++j){
        for(int k = 0; j < lim[2]; ++k){
          double e = ( vibf[0] * (i + 0.5) + vibf[1] * (j + 0.5) + vibf[2] * (k + 0.5) ) * h * c;
          //std::cout << i << " " << j << " " << k << "\n";
          if(e < de) { vibl.push_back(e); }
          else { break; }
        }
      }
    }
    double Z_vib = 0;
    for(const double& lvl : vibl){
      Z_vib += std::exp(-lvl / (K * T));
    }
    std::cout << Z_vib;
    return Z_vib;
  }
  double Z_vib = 1;
  for(const auto& freq : comp.vibrational_freq){
    double denom = 1 - std::exp(-h * c * freq / (K * T));
    Z_vib *= (denom > 0) ? 1/denom : 1.0;
  }
  return Z_vib;
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
