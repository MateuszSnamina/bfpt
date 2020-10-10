#include<monostar_app/interpret_orbital_theta_string.hpp>

#include<monostar_app/interpret_number_string.hpp>

namespace monostar_app {

utility::Result<std::optional<double>, std::domain_error> interpret_orbital_theta_string(const std::string& s) {
    using ResultT = utility::Result<std::optional<double>, std::domain_error>;
    if (s == "auto") {
        return ResultT::Ok(std::nullopt);
    } else {
        if (const auto& _ = interpret_double_string(s)) {
            return ResultT::Ok(_.unwrap());
        } else {
            std::string message{"Orbital theta string should be a number or `auto`."};
            std::domain_error error{message};
            return ResultT::Err(error);
        }
    }
}

} // end of namespace monostar_app
