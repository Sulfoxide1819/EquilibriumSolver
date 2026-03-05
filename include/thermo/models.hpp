#pragma once
class VibrationalModel {
public:
  virtual ~VibrationalModel() = default;
  virtual double compute(double T) const = 0;
};

class HarmonicOscillator: public VibrationalModel {
  std::vector<double> freqs_;
public:
  HarmonicOscillator(const std::vector<double>& freqs_): freqs_(freqs_){}
  double compute(double T) const override;
};

class AnharmonicOscillator: public VibrationalModel {
public:
  AnharmonicOscillator(double diss_energy, const std::vector<double>& freqs_, std::vector<double> omega_x):
    freqs_(freqs_),
    diss_energy(diss_energy),
    omega_x(omega_x) {}
  double compute(double T) const override;
private:
  double diss_energy;
  std::vector<double> freqs_;
  std::vector<double> omega_x;
};

class CO2_model : public VibrationalModel {
public:
  CO2_model(double diss_energy, const std::vector<double>& freqs_, const Eigen::Matrix3d& anh_consts2) :
    freqs_(freqs_),
    aharmonic_constants2(anh_consts2),
    diss_energy(diss_energy) {}

  double compute(double T) const override;
private:
  double diss_energy;//per molecule
  std::vector<double> freqs_;
  Eigen::Matrix3d aharmonic_constants2;//triangular matrix
};
