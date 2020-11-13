#pragma once

#include <utility/result.hpp>

#include <stdexcept>
#include <optional>
#include <string>

namespace monostar_app {

utility::Result<std::optional<double>, std::domain_error> interpret_orbital_theta_string(const std::string&);

}  // end of namespace monostar_app
