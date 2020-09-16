#pragma once

#include <model_monostar/enum_model_type.hpp>

#include <utility/result.hpp>

#include <map>

#include <stdexcept>
#include <string>

extern const std::map<std::string, ModelType> interpret_model_type_string_map;

utility::Result<ModelType, std::domain_error> interpret_model_type_string(const std::string&);
