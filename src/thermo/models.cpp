#include <stdexcept>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <memory>
#include "core/constants.hpp"
#include "thermo/models.hpp"

#include <iostream>

double HarmonicOscillator::compute(double T) const {
  double Z_vib = 1;
  for(const double& freq : freqs_){
    double denom = 1 - std::exp(-h * c * freq / (K * T));
    Z_vib *= (denom > 0) ? 1.0/denom : 1.0;
  }
  return Z_vib;
}

double AnharmonicOscillator::compute(double T) const {
    if(freqs_.size() > 1) throw std::runtime_error("AnharmonicOscillator model do not support non-diatomic molecule");
    size_t v = 0;
    double e_base = freqs_[0] * 0.5 - omega_x[0] * 0.25;
    double e = 0;
    double Z_vib = 0;
    size_t v_max = 39;

    while(true){
      e = freqs_[0] * (v + 0.5) - omega_x[0] * (v + 0.5) * (v + 0.5);
      e = (e - e_base) * h * c;
      if(e < diss_energy) {
        double e_level = std::exp(-e / (K * T));
        if (e_level < 1e-6 || v > v_max) { break;}
        Z_vib += e_level;
        ++v; 
      } else { break; }
    }
    return Z_vib;
  }
  

double CO2_model::compute(double T) const {
  std::vector<unsigned> lim = {34, 67, 20};
  double cut_energy = 8.83859e-19;
  double Z_vib = 0;

  for(int i = 0; i < lim[0]; ++i){
    for(int j = 0; j < lim[1]; ++j){
      for(int k = 0; k < lim[2]; ++k){
        Eigen::Vector3d l;//a vector of quantum number
        l << i, j, k;

        double e_har = freqs_[0] * i + freqs_[1] * j + freqs_[2] * k;//harmonic term
        double e_anh = l.transpose() * aharmonic_constants2 * l;//anharmonic term
        double e = (e_har + e_anh) * h * c;

        if(e < cut_energy) { Z_vib += (j + 1) * std::exp(-e / (K * T)); }
        else { break; }
      }
    }
  }
  return Z_vib;
} 
