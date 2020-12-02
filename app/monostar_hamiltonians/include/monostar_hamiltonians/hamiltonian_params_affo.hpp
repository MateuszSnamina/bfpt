#pragma once

#include <monostar_hamiltonians/hamiltonian_params_fo.hpp>
#include <monostar_hamiltonians/hamiltonian_params_af_fm.hpp>

#include <set>
#include <string>

// #######################################################################
// ## HamiltonianParamsAffo                                             ##
// #######################################################################

namespace monostar_hamiltonians {

class HamiltonianParamsAffo {
   public:
    class Builder {
       public:
        Builder set_s_coef(double);
        Builder set_ss_coef(double);
        Builder set_tau_z_coef(double);
        Builder set_tau_minus_coef(double);
        Builder set_Pzz_coef(double);
        Builder set_Pxz_coef(double);
        Builder set_Pxx_coef(double);
        Builder set_ss_Pzz_coef(double);
        Builder set_ss_Pxz_coef(double);
        Builder set_ss_Pxx_coef(double);
        HamiltonianParamsAffo build() const;

       private:
        double _s_coef = 0.0;
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
    friend HamiltonianParamsAffo Builder::build() const;
    double get_s_coef() const;
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
    double get_site_energy(double theta, double average_ss) const;
    double get_site_energy_derivative(double theta, double average_ss) const;
    double get_site_energy_derivative2(double theta, double average_ss) const;
    double get_site_energy_derivative3(double theta, double average_ss) const;
    double get_site_energy_derivative4(double theta, double average_ss) const;
    monostar_hamiltonians::HamiltonianParamsFo average_out_spins_1(double average_s) const;
    monostar_hamiltonians::HamiltonianParamsFo average_out_spins_12(double average_s, double average_ss) const;
    monostar_hamiltonians::HamiltonianParamsAfFm average_out_orbitals_1(double theta) const;
    monostar_hamiltonians::HamiltonianParamsAfFm average_out_orbitals_12(double average_tau_minus, double average_tau_z, double average_Pzz, double average_Pzx_sum_P_xz, double average_Pxx) const;
    std::set<double> get_theta_opt(double average_ss) const;
    std::set<double> get_theta_opt_numerical(double average_ss) const;
    utility::Result<std::set<double>, monostar_hamiltonians::NoKnownAnalyticalSolutionError> get_theta_opt_analytical(double average_ss) const;
    std::string string_repr_in_orbital_operators() const;
    std::string string_repr_in_trigonometric_functions() const;

   private:
    HamiltonianParamsAffo(
        double s_coef, double ss_coef,
        double tau_z_coef, double tau_minus_coef,
        double Pzz_coef, double Pxz_coef, double Pxx_coef,
        double ss_Pzz_coef, double ss_Pxz_coef, double ss_Pxx_coef);
    double _s_coef = 0.0;
    double _ss_coef = 0.0;
    monostar_hamiltonians::HamiltonianParamsFo _hamiltonian_params_fo;
    monostar_hamiltonians::HamiltonianParamsFo _hamiltonian_params_ss_fo;
};

}  // end of namespace monostar_hamiltonians
