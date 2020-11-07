#pragma once

#include <so_app/enum_model_type.hpp>

#include <utility/result.hpp>

#include <map>
#include <stdexcept>
#include <string>

namespace so_app {

extern const std::map<std::string, ModelType> interpret_model_type_string_map;

utility::Result<ModelType, std::domain_error> interpret_model_type_string(const std::string&);

}  // end of namespace so_app
