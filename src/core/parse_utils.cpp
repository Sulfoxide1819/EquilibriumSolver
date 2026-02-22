#include <string>
#include <nlohmann/json.hpp>
#include <vector>
#include "core/types.hpp"
#include "core/parse_utils.hpp"
#include <stdexcept>
#include <set>

using json = nlohmann::json;


  void MixtureBuild::pull_properties(const json& data, Component& comp){
    if(!data.contains(comp.name)) throw std::runtime_error("No data for " + comp.name + "!");

    auto& comp_data = data[comp.name];

    comp.molar_mass = comp_data.at("molar_mass").get<double>();
    comp.is_atomic = comp_data.at("is_atomic").get<bool>();
    if(!comp.energy_levels.empty()) {
      for(const auto& lvl : comp_data["energy_levels"]){
        Component::EnergyLevel level;
        level.state = lvl.at("state").get<std::string>();
        level.energy = lvl.at("T_i").get<double>();
        level.degeneracy = lvl.at("p_i").get<int>();
        comp.energy_levels.push_back(level);
      } 
     } 
     if(!comp.is_atomic) {
          comp.dissociation_energy = comp_data["dissociation_energy"].get<double>();
          comp.vibrational_freq = comp_data["vibrational_freq"].get<std::vector<double>>();
          comp.rotational_const = comp_data["rotational_const"].get<double>();
          comp.symmetry_factor = comp_data["symmetry_factor"].get<int>();
          if (comp_data.contains("anharmonic_constants_2")){
            size_t size = comp.vibrational_freq.size();
            Eigen::MatrixXd matr= Eigen::MatrixXd::Zero(size, size);
            size_t it = 0;
            for(int i = 0; i < size; ++i){
              for(int j = i; j < size; ++j){
                matr(i, j) = comp_data["anharmonic_constants_2"].at(it).get<double>();
                ++it;
                comp.anharmonic_constants_2 = matr;
              }

            }

          }
     }
    }
  
  std::vector<Element> MixtureBuild::get_elements(const json& data, const std::vector<Component>& comps) {
    std::set<std::string> set_el; //set of elements
    for(const auto& comp : comps){
      if(!data.at(comp.name).contains("stoichiometry")) throw std::runtime_error("No stoichiometry for " + comp.name);
      for(const auto& el : data[comp.name]["stoichiometry"].items()) {
        if(!set_el.count(el.key())) { set_el.insert(el.key()); }
      }
    }
    std::vector<Element> elements;
    for(const auto& el : set_el){
      elements.push_back(Element(el));
    }
    return elements;
  }
  

  void MixtureBuild::get_stoichiometry(const json& data, Element& el, const  std::vector<Component>& comps){
    //el.stoichiometry.resize(comp.size());
    for(const auto& comp : comps){
  
      int value = 0;
      if(data[comp.name]["stoichiometry"].contains(el.name)) { value = data[comp.name]["stoichiometry"].at(el.name).get<int>(); }
      el.stoichiometry.push_back(value);
    }
  }

  int MixtureBuild::pull_thermo(const json& data, Component& comp, double T) {
    if(data.contains(comp.name)){
      const auto& comp_data = data[comp.name];
      size_t id = 0;
      for(const auto& temp : comp_data["T"]){
        if(T > temp){ ++id; continue; }
        --id; break;
      }
      comp.reduced_gibbs_energy = comp_data["Phi"].at(id).get<double>();
      comp.d_entalpy_0 = comp_data["H"].at(id).get<double>();
      return 0;
    }
    return 1;

  }

 

//MixtureBuild

