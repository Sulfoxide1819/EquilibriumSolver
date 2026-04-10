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
  return K * T / (c * h * comp.symmetry_factor * comp.rotational_const);
}

double StatSum::vibrational(const Component& comp, double T) {
  if(comp.is_atomic) return 1.0;
  double Z_vib = comp.vib_model->compute(T);
  return Z_vib;
}
/*double StatSum::vibrational(const Component& comp, double T) {
  if (comp.is_atomic) return 1.0;
  const std::vector<double>& vibf = comp.vibrational_freq;
  if(vibf.size() > 1) {
    std::vector<double> vibl; //energy of vibrational levels
    double de = comp.dissociation_energy / NA;;
    std::vector<int> lim = {34, 67, 20};
   // std::vector<size_t> n(0,vibf.size());
    //size_t index = vibf.size() - 1;
    /-/while(index >= 0){
      double e = 0;
      double e_har = ( vibf[0] * (i + 0.5) + vibf[1] * (j + 1.0) + vibf[2] * (k + 0.5) ) * h * c; //harmonic_term
      double e_anh = 0;
      const auto& matr = comp.anharmonic_constants_2;
      Eigen::VectorXd l(3);
      l << i + 0.5, j + 1.0, k + 0.5;
      e_anh = (l.transpose() * matr * l);
      e = e_har + e_anh * h * c;
      if(e < de){
        vibl.push_back(e);
        size_t size = vibf.size();
        ++n[size - 1];
        index = size - 1;
      } else {
        n[index] = 0;
        --index;
        ++n[index];
      }
    }/-/
    
    double e_base_har = ( vibf[0] * 0.5 + vibf[1] * 1.0 + vibf[2] * 0.5 ) * h * c;
    const auto& matr = comp.anharmonic_constants_2;
    double e_base_anh = (Eigen::VectorXd(3) << 0.5, 1.0, 0.5).finished().transpose() * matr * (Eigen::VectorXd(3) << 0.5, 1.0, 0.5).finished();
    std::cout << matr << "\n";
    double Z_vib = 0;
    for(int i = 0; i < lim[0]; ++i){
      for(int j = 0; j < lim[1]; ++j){
        for(int k = 0; k < lim[2]; ++k){
          double e = 0;
          double e_har = ( vibf[0] * i + vibf[1] * j + vibf[2] * k  ) * h * c; //harmonic_term
          double e_anh = 0;
          Eigen::VectorXd l(3);
          l << i, j, k;
          e_anh = (l.transpose() * matr * l);
          e = e_har + e_anh * h * c;
          if(e < de) { //vibl.push_back(e); 
            Z_vib += (j + 1) * std::exp(-e/(K * T));
          }
          else { break; }
        }
      }
    }
    double Z_vib = 0;
    for(const double& lvl : vibl){
      Z_vib += std::exp(-lvl / (K * T));
    }
    std::cout << "vib_Z (CO2):" << Z_vib << "\n";
    return Z_vib;
    //return 43.0131;
  }
 // double Z_vib = 0;
  if(!comp.anharmonic_constants_2.empty()) { 
    double de = comp.dissociation_energy / NA;
    double omega = comp.vibrational_freq[0];
    double omega_x = comp.anharmonic_constants_2(0,0);
    size_t v = 0;
    do {
      double e = h * c * (omega * (v + 0.5) - omega_x * (v + 0.5) * (v + 0.5));
      ++v;
      Z_vib += std::exp(-e/ (K * T));
    } while(e < de)
    return Z_vib;

  }
  double Z_vib = 1;
  for(const auto& freq : comp.vibrational_freq){
    double denom = 1 - std::exp(-h * c * freq / (K * T));
    Z_vib *= (denom > 0) ? 1/denom : 1.0;
  }
  std::cout << "vib_Z: " << Z_vib << "\n";
  return Z_vib;
}*/

/*double vibrational_adv(const Component& comp, const std::vector<double>& freq, double T, size_t lay){
  
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
/*  if (comp.name == "CO2") {
    double res = std::exp(fromJANAF(comp, T, 0));
    std::cout << "Z total CO2:" << res << "\n";
    return res;}
  if (comp.name == "CO") {
    double res = std::exp(fromJANAF(comp, T, 1));
    return res;}
  if (comp.name == "C") {
    double res = std::exp(fromJANAF(comp, T, 2));
    return res;}
  if (comp.name == "O") {
    double res = std::exp(fromJANAF(comp, T, 3));
    return res;}
  if (comp.name == "O2") {
    double res = std::exp(fromJANAF(comp, T, 4));
    return res;}*/
  double exp_term = std::exp(comp.dissociation_energy / (K * T * NA));
  double Z_total = translational(comp, T) * rotational(comp, T) * vibrational(comp, T) * electronic(comp, T) * exp_term;
  //std::cout << "Z_total(" + comp.name + "): " << Z_total << "\n";
  return Z_total;
}

double StatSum::fromReducedGibbs(const Component& comp, double T){
  const double& Phi = comp.reduced_gibbs_energy;
  const double& dH = comp.dissociation_energy;
  double ln_Z = Phi / (NA * K) + std::log(NA * p_0 / (NA * K * T)) + dH / (NA * K * T);
  return ln_Z;

}

double fromJANAF(const Component& comp, double T, int n){
  std::vector<double> S_0 = {334.169, 273.623, 206.322, 209.704, 284.466};
  std::vector<double> dH = {152.852, 93.546, 56.689, 56.574, 98.013};
  double d_fH_0 = comp.dissociation_energy;
  constexpr double R = K * NA;
  double lnZ = S_0[n] / R - dH[n] * 1e3 / (R * T) + d_fH_0 / (R * T) + std::log(NA * p_0 / (R * T));
  return lnZ;
}



