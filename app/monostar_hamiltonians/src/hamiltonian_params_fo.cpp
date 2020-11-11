#include <monostar_hamiltonians/hamiltonian_params_fo.hpp>
#include <monostar_hamiltonians/hamiltonian_params_fo_helpers.hpp>
#include<monostar_hamiltonians/hamiltonian_params_fo_to_classic_energy_function.hpp>

#include <sstream>
#include <iomanip>

// #######################################################################
// ## HamiltonianParamsFo                                               ##
// #######################################################################

namespace monostar_hamiltonians {

HamiltonianParamsFo::Builder HamiltonianParamsFo::Builder::set_tau_z_coef(double tau_z_coef) {
    _tau_z_coef = tau_z_coef;
    return *this;
}

HamiltonianParamsFo::Builder HamiltonianParamsFo::Builder::set_tau_minus_coef(double tau_minus_coef) {
    _tau_munis_coef = tau_minus_coef;
    return *this;
}

HamiltonianParamsFo::Builder HamiltonianParamsFo::Builder::set_Pzz_coef(double Pzz_coef) {
    _Pzz_coef = Pzz_coef;
    return *this;
}

HamiltonianParamsFo::Builder HamiltonianParamsFo::Builder::set_Pxz_coef(double Pxz_coef) {
    _Pxz_coef = Pxz_coef;
    return *this;
}

HamiltonianParamsFo::Builder HamiltonianParamsFo::Builder::set_Pxx_coef(double Pxx_coef) {
    _Pxx_coef = Pxx_coef;
    return *this;
}

HamiltonianParamsFo::Builder HamiltonianParamsFo::Builder::set_free_coef(double free_coef) {
    _free_coef = free_coef;
    return *this;
}

HamiltonianParamsFo HamiltonianParamsFo::Builder::build() const {
    return HamiltonianParamsFo(_tau_z_coef, _tau_munis_coef, _Pzz_coef, _Pxz_coef, _Pxx_coef, _free_coef);
}

// #######################################################################

HamiltonianParamsFo::HamiltonianParamsFo(double tau_z_coef, double tau_minus_coef, double Pzz_coef, double Pxz_coef, double Pxx_coef, double free_coef)
    : _tau_z_coef(tau_z_coef),
      _tau_minus_coef(tau_minus_coef),
      _Pzz_coef(Pzz_coef),
      _Pxz_coef(Pxz_coef),
      _Pxx_coef(Pxx_coef),
      _free_coef(free_coef){
}

double HamiltonianParamsFo::get_tau_z_coef() const {
    return _tau_z_coef;
}

double HamiltonianParamsFo::get_tau_minus_coef() const {
    return _tau_minus_coef;
}

double HamiltonianParamsFo::get_Pzz_coef() const {
    return _Pzz_coef;
}

double HamiltonianParamsFo::get_Pxz_coef() const {
    return _Pxz_coef;
}

double HamiltonianParamsFo::get_Pxx_coef() const {
    return _Pxx_coef;
}

double HamiltonianParamsFo::get_free_coef() const {
    return _free_coef;
}

double HamiltonianParamsFo::get_site_energy(double theta) const {
    return hamiltonian_params_fo_to_classic_energy_function(*this).get_value(theta);
}

double HamiltonianParamsFo::get_site_energy_derivative(double theta) const {
    return hamiltonian_params_fo_to_classic_energy_function(*this).get_derivative_value(theta);
}

double HamiltonianParamsFo::get_site_energy_derivative2(double theta) const {
    return hamiltonian_params_fo_to_classic_energy_function(*this).get_derivative2_value(theta);
}

double HamiltonianParamsFo::get_site_energy_derivative3(double theta) const {
    return hamiltonian_params_fo_to_classic_energy_function(*this).get_derivative3_value(theta);
}

double HamiltonianParamsFo::get_site_energy_derivative4(double theta) const {
    return hamiltonian_params_fo_to_classic_energy_function(*this).get_derivative4_value(theta);
}

std::set<double> HamiltonianParamsFo::get_theta_opt() const {
    return hamiltonian_params_fo_to_classic_energy_function(*this).get_minimum_argument();
}

std::set<double> HamiltonianParamsFo::get_theta_opt_numerical() const {
    return hamiltonian_params_fo_to_classic_energy_function(*this).get_minimum_argument_numerical();
}

utility::Result<std::set<double>, NoKnownAnalyticalSolutionError> HamiltonianParamsFo::get_theta_opt_analytical() const {
    return hamiltonian_params_fo_to_classic_energy_function(*this).get_minimum_argument_analytical();
}

std::string HamiltonianParamsFo::string_repr_in_orbital_operators() const {
    std::ostringstream ss;
    ss << std::showpos;
    if (_tau_z_coef != 0) {
        ss << _tau_z_coef << "·Σᵢ[τᶻᵢ]";
    }
    if (_tau_minus_coef != 0) {
        ss << _tau_minus_coef << "·Σᵢ[τ⁻ᵢ]";
    }
    if (_Pzz_coef != 0) {
        ss << _Pzz_coef << "·Σᵢⱼ[PᶻᵢPᶻⱼ]";
    }
    if (_Pxz_coef != 0) {
        ss << _Pxz_coef << "·Σᵢⱼ[PˣᵢPᶻⱼ + PᶻᵢPˣⱼ]";
    }
    if (_Pxx_coef != 0) {
        ss << _Pxx_coef << "·Σᵢⱼ[PˣᵢPˣⱼ]";
    }
    if (_free_coef != 0) {
        ss << _free_coef;
    }
    if (_tau_z_coef == 0 && _tau_minus_coef == 0 && _Pzz_coef == 0 && _Pxz_coef == 0 && _Pxx_coef == 0 && _free_coef == 0) {
        ss << 0.0;
    }
    return ss.str();
}

std::string HamiltonianParamsFo::string_repr_in_trigonometric_functions() const {
    const auto classic_energy_function = hamiltonian_params_fo_to_classic_energy_function(*this);
    const double sin_coef = classic_energy_function.get_sin_coef();
    const double cos_coef = classic_energy_function.get_cos_coef();
    const double sqcos_coef = classic_energy_function.get_sqcos_coef();
    const double free_coef = classic_energy_function.get_free_coef();
    std::ostringstream ss;
    ss << std::showpos;
    if (sin_coef != 0) {
        ss << sin_coef << "·sin(θ)";
    }
    if (cos_coef != 0) {
        ss << cos_coef << "·cos(θ)";
    }
    if (sqcos_coef != 0) {
        ss << sqcos_coef << "·cos²(θ)";
    }
    if (free_coef != 0) {
        ss << free_coef;
    }
    if (sin_coef == 0 && cos_coef == 0 && sqcos_coef == 0 && free_coef == 0) {
        ss << 0.0;
    }
    return ss.str();
}

}  // end of namespace monostar_hamiltonians
