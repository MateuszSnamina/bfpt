#pragma once

#include <utility/result.hpp>

#include <stdexcept>
#include <optional>
#include <string>

namespace so_app {

utility::Result<std::optional<unsigned>, std::domain_error> interpret_n_max_site_excitations_string(const std::string&);

}  // end of namespace so_app
