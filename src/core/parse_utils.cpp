#include <string>
#include <nlohmann/json.hpp>
#include <vector>
#include "core/types.hpp"
#include "core/parse_utils.hpp"
#include <stdexcept>

using json = nlohmann::json;


  void MixtureBuild::pull_properties(const json& data, Component& comp){
    if(!data.contains(comp.name)) throw std::runtime_error("No data for " + comp.name + "!");

    auto& comp_data = data[comp.name];

    comp.molar_mass = comp_data.at("molar_mass").get<double>();
    comp.is_atomic = comp_data.at("is_atomic").get<bool>();
    if(comp.is_atomic) {
      for(const auto& lvl : comp_data["energy_levels"]){
        Component::EnergyLevel level;
        level.state = lvl.at("state").get<std::string>();
        level.energy = lvl.at("T_i").get<double>();
        level.degeneracy = lvl.at("p_i").get<int>();
        comp.energy_levels.push_back(level);
      } 
     } else {
          comp.dissociation_energy = comp_data["dissociation_energy"].get<double>();
          comp.vibrational_freq = comp_data["vibrational_freq"].get<std::vector<double>>();
          comp.rotational_const = comp_data["rotational_const"].get<double>();
          comp.symmetry_factor = comp_data["symmetry_factor"].get<int>();
     }
    }
  

  

  void MixtureBuild::get_stoichiometry(const json& data, Element& el, const  std::vector<Component>& comps){
    //el.stoichiometry.resize(comp.size());
    for(const auto& comp : comps){
  
      int value = 0;
      if(data[comp.name]["stoichiometry"].contains(el.name)) { value = data[comp.name]["stoichiometry"].at(el.name).get<int>(); }
      el.stoichiometry.push_back(value);
    }


  }

//MixtureBuild

