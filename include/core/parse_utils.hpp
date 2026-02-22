#pragma once

#include "core/types.hpp"
#include <nlohmann/json.hpp>

namespace MixtureBuild {

void pull_properties(const nlohmann::json& data, Component& comp);

std::vector<Element> get_elements(const nlohmann::json& data, const std::vector<Component>& comps);
void get_stoichiometry(const nlohmann::json& data, Element& el, const std::vector<Component>& comps);

int pull_thermo(const nlohmann::json& data, Component& comps, double T);
}
