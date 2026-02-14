#pragma once

#include "core/types.hpp"
#include <nlohmann/json.hpp>

namespace MixtureBuild {

void pull_properties(const nlohmann::json& data, Component& comp);

void get_stoichiometry(const nlohmann::json& data, Element& el, const std::vector<Component>& comps);

}
