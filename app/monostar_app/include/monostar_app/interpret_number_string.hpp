#pragma once

// UTILITY
#include <utility/result.hpp>
// STD:
#include <string>
#include <optional>
#include <stdexcept>

//**********************************************************************
//***   interpret_{int,unsigned,double}_string                       ***
//**********************************************************************

namespace monostar_app {

utility::Result<double, std::domain_error> interpret_double_string(const std::string&);
utility::Result<int, std::domain_error> interpret_int_string(const std::string&);
utility::Result<unsigned, std::domain_error> interpret_unsigned_string(const std::string&);

}  // namespace monostar_app
