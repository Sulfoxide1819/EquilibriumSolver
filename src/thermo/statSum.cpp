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
  std::cout << "Z_tr ( " << comp.name << " ): " << Z_tr << "\n";
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
        for(int k = 0; k < lim[2]; ++k){
          double e = 0;
          double e_har = ( vibf[0] * (i + 0.5) + vibf[1] * (j + 1.0) + vibf[2] * (k + 0.5) ) * h * c; //harmonic_term
          double e_anh = 0;
          const auto& matr = comp.anharmonic_constants_2;
          Eigen::VectorXd l(3);
          l << i + 0.5, j + 1.0, k + 0.5;
          e_anh = (l.transpose() * matr * l);
          e = e_har + e_anh * h * c;
          if(e < de) { vibl.push_back(e); }
          else { break; }
        }
      }
    }
    double Z_vib = 0;
    for(const double& lvl : vibl){
      Z_vib += std::exp(-lvl / (K * T));
    }
    std::cout << "vib_Z: " << Z_vib << "\n";
    //return Z_vib;
    return 43.0131;
  }
  double Z_vib = 1;
  for(const auto& freq : comp.vibrational_freq){
    double denom = 1 - std::exp(-h * c * freq / (K * T));
    Z_vib *= (denom > 0) ? 1/denom : 1.0;
  }
  std::cout << "vib_Z: " << Z_vib << "\n";
  return Z_vib;
}

/*double vibrational_adv(const Component& comp, double T, size_t m){
  if(m > 0) {
    for(int i = 0; ;++i){
      double e = vibrational_adv(comp, T, m - 1);
    }
  else {
    
  }
}*/

double StatSum::electronic(const Component& comp, double T){
  if(comp.energy_levels.empty()){ return 1.0; }
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

double StatSum::fromReducedGibbs(const Component& comp, double T){
  const double& Phi = comp.reduced_gibbs_energy;
  const double& dH = comp.dissociation_energy;
  double ln_Z = Phi / (NA * K) + std::log(NA * p_0 / (NA * K * T)) + dH / (NA * K * T);
  return ln_Z;

}





