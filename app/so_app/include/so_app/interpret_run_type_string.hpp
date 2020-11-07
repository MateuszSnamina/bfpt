#pragma once

#include <so_app/enum_run_type.hpp>

#include <utility/result.hpp>

#include <map>

#include <stdexcept>
#include <string>

namespace so_app {

extern const std::map<std::string, RunType> interpret_run_type_string_map;

utility::Result<RunType, std::domain_error> interpret_run_type_string(const std::string&);

}  // end of namespace so_app
