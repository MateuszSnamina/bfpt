#include <so_app/interpret_n_max_site_excitations_string.hpp>

#include <so_app/interpret_number_string.hpp>

namespace so_app {

utility::Result<std::optional<double>, std::domain_error> interpret_n_max_site_excitations_string(const std::string& s) {
    using ResultT = utility::Result<std::optional<double>, std::domain_error>;
    if (s == "nolimit") {
        return ResultT::Ok(std::nullopt);
    } else {
        if (const auto& _ = interpret_double_string(s)) {
            return ResultT::Ok(_.unwrap());
        } else {
            std::string message{"Max site excitations string should be a number or `nolimit`."};
            std::domain_error error{message};
            return ResultT::Err(error);
        }
    }
}

}  // end of namespace so_app
