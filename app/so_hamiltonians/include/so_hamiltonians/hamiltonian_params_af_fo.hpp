#pragma once

#include <monostar_hamiltonians/hamiltonian_params_fo.hpp>

#include <set>
#include <string>

// #######################################################################
// ## HamiltonianParamsFo                                               ##
// #######################################################################

namespace so_hamiltonians {

class HamiltonianParamsAfFo {
   public:
    class Builder {
       public:
        Builder set_ss_coef(double);
        Builder set_tau_z_coef(double);
        Builder set_tau_minus_coef(double);
        Builder set_Pzz_coef(double);
        Builder set_Pxz_coef(double);
        Builder set_Pxx_coef(double);
        Builder set_ss_Pzz_coef(double);
        Builder set_ss_Pxz_coef(double);
        Builder set_ss_Pxx_coef(double);
        HamiltonianParamsAfFo build() const;

       private:
        double _ss_coef = 0.0;
        double _tau_z_coef = 0.0;
        double _tau_munis_coef = 0.0;
        double _Pzz_coef = 0.0;
        double _Pxz_coef = 0.0;
        double _Pxx_coef = 0.0;
        double _ss_Pzz_coef = 0.0;
        double _ss_Pxz_coef = 0.0;
        double _ss_Pxx_coef = 0.0;
    };
    friend HamiltonianParamsAfFo Builder::build() const;
    HamiltonianParamsAfFo() = default;
    double get_ss_coef() const;
    double get_tau_z_coef() const;
    double get_tau_minus_coef() const;
    double get_Pzz_coef() const;
    double get_Pxz_coef() const;
    double get_Pxx_coef() const;
    double get_ss_Pzz_coef() const;
    double get_ss_Pxz_coef() const;
    double get_ss_Pxx_coef() const;
    monostar_hamiltonians::HamiltonianParamsFo get_hamiltonian_params_fo() const;
    monostar_hamiltonians::HamiltonianParamsFo get_hamiltonian_params_ss_fo() const;
    double get_site_energy(double theta) const;
    double get_site_energy_derivative(double theta) const;
    double get_site_energy_derivative2(double theta) const;
    double get_site_energy_derivative3(double theta) const;
    double get_site_energy_derivative4(double theta) const;
    std::set<double> get_theta_opt() const;
    std::set<double> get_theta_opt_numerical() const;
    utility::Result<std::set<double>, monostar_hamiltonians::NoKnownAnalyticalSolutionError> get_theta_opt_analytical() const;
    std::string string_repr_in_orbital_operators() const;
    std::string string_repr_in_trigonometric_functions() const;

   private:
    HamiltonianParamsAfFo(
            double ss_coef,
            double tau_z_coef, double tau_minus_coef, double Pzz_coef, double Pxz_coef, double Pxx_coef,
            double ss_Pzz_coef, double ss_Pxz_coef, double ss_Pxx_coef);
    double _ss_coef = 0.0;
    monostar_hamiltonians::HamiltonianParamsFo _hamiltonian_params_fo;
    monostar_hamiltonians::HamiltonianParamsFo _hamiltonian_params_ss_fo;
};

}  // end of namespace so_hamiltonians
