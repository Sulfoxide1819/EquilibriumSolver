#pragma once
#include "core/types.hpp"
class StatSum {
public:
  static double translational(const Component& comp, double T);
  static double rotational(const Component& comp, double T);
  static double vibrational(const Component& comp, double T);
  static double electronic(const Component& comp, double T);
  static double total(const Component& comp, double T);
};
