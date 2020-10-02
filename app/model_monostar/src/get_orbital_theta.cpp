#include<model_monostar/get_orbital_theta.hpp>

#include <boost/range/adaptor/filtered.hpp>

#include <string>

double get_orbital_theta(
        const HamiltonianFoParams& hamiltonian_fo_params,
        std::optional<double> user_defined_overrule) {
    if (user_defined_overrule) {
        return *user_defined_overrule;
    } else {
        const auto theta_opt_set = hamiltonian_fo_params.get_theta_opt();
        const auto positive_theta_opt_set = theta_opt_set | boost::adaptors::filtered([](double theta){return theta >=0;});
        if (!boost::empty(positive_theta_opt_set)) {
            return *positive_theta_opt_set.begin();
        } else if (!boost::empty(theta_opt_set)) {
            return *theta_opt_set.crbegin();
        } else {
            const std::string message{"Cannot select orbital theta as for the given hamiltonian params there is no optinal theta."};
            const std::runtime_error error{message};
            throw error;
        }
    }
}
