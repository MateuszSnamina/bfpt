#include <so_hamiltonians/hamiltonian_params_af_fo.hpp>

#include <monostar_hamiltonians/hamiltonian_params_fo_helpers.hpp>
#include <monostar_hamiltonians/hamiltonian_params_fo_to_classic_energy_function.hpp>

#include <sstream>
#include <iomanip>

// #######################################################################
// ## hamiltonian_params_af_fo_to_classic_energy_function               ##
// #######################################################################

namespace so_hamiltonians {

monostar_hamiltonians::AcosPlusBsinPlusCsqcosPlusZ hamiltonian_params_af_fo_to_classic_energy_function(HamiltonianParamsAfFo params) {
    double SS_average = -4.0; //TODO: make SS_average an arg.
    const monostar_hamiltonians::AcosPlusBsinPlusCsqcosPlusZ fun1 =
            monostar_hamiltonians::hamiltonian_params_fo_to_classic_energy_function(params.get_hamiltonian_params_fo());
    const monostar_hamiltonians::AcosPlusBsinPlusCsqcosPlusZ fun2 =
            monostar_hamiltonians::hamiltonian_params_fo_to_classic_energy_function(params.get_hamiltonian_params_ss_fo());
    return fun1 + (SS_average * fun2);//TODO finish
}

}  // end of namespace so_hamiltonians

// #######################################################################
// ## HamiltonianParamsAfFo                                             ##
// #######################################################################

namespace so_hamiltonians {

HamiltonianParamsAfFo::Builder HamiltonianParamsAfFo::Builder::set_ss_coef(double ss_coef) {
    _ss_coef = ss_coef;
    return *this;
}

HamiltonianParamsAfFo::Builder HamiltonianParamsAfFo::Builder::set_tau_z_coef(double tau_z_coef) {
    _tau_z_coef = tau_z_coef;
    return *this;
}

HamiltonianParamsAfFo::Builder HamiltonianParamsAfFo::Builder::set_tau_minus_coef(double tau_minus_coef) {
    _tau_munis_coef = tau_minus_coef;
    return *this;
}

HamiltonianParamsAfFo::Builder HamiltonianParamsAfFo::Builder::set_Pzz_coef(double Pzz_coef) {
    _Pzz_coef = Pzz_coef;
    return *this;
}

HamiltonianParamsAfFo::Builder HamiltonianParamsAfFo::Builder::set_Pxz_coef(double Pxz_coef) {
    _Pxz_coef = Pxz_coef;
    return *this;
}

HamiltonianParamsAfFo::Builder HamiltonianParamsAfFo::Builder::set_Pxx_coef(double Pxx_coef) {
    _Pxx_coef = Pxx_coef;
    return *this;
}

HamiltonianParamsAfFo::Builder HamiltonianParamsAfFo::Builder::set_ss_Pzz_coef(double ss_Pzz_coef) {
    _ss_Pzz_coef = ss_Pzz_coef;
    return *this;
}

HamiltonianParamsAfFo::Builder HamiltonianParamsAfFo::Builder::set_ss_Pxz_coef(double ss_Pxz_coef) {
    _ss_Pxz_coef = ss_Pxz_coef;
    return *this;
}

HamiltonianParamsAfFo::Builder HamiltonianParamsAfFo::Builder::set_ss_Pxx_coef(double ss_Pxx_coef) {
    _ss_Pxx_coef = ss_Pxx_coef;
    return *this;
}

HamiltonianParamsAfFo HamiltonianParamsAfFo::Builder::build() const {
    return HamiltonianParamsAfFo(
                _ss_coef,
                _tau_z_coef, _tau_munis_coef, _Pzz_coef, _Pxz_coef, _Pxx_coef,
                _ss_Pzz_coef, _ss_Pxz_coef, _ss_Pxx_coef);
}

// #######################################################################

HamiltonianParamsAfFo::HamiltonianParamsAfFo(
        double ss_coef,
        double tau_z_coef, double tau_minus_coef, double Pzz_coef, double Pxz_coef, double Pxx_coef,
        double ss_Pzz_coef, double ss_Pxz_coef, double ss_Pxx_coef)
    : _ss_coef(ss_coef),
      _hamiltonian_params_fo(monostar_hamiltonians::HamiltonianParamsFo::Builder()
                            .set_tau_z_coef(tau_z_coef)
                            .set_tau_minus_coef(tau_minus_coef)
                            .set_Pzz_coef(Pzz_coef)
                            .set_Pxz_coef(Pxz_coef)
                            .set_Pxx_coef(Pxx_coef)
                            .build()),
      _hamiltonian_params_ss_fo(monostar_hamiltonians::HamiltonianParamsFo::Builder()
                               .set_tau_z_coef(0.0)
                               .set_tau_minus_coef(0.0)
                               .set_Pzz_coef(ss_Pzz_coef)
                               .set_Pxz_coef(ss_Pxz_coef)
                               .set_Pxx_coef(ss_Pxx_coef)
                               .build()) {
}

double HamiltonianParamsAfFo::get_ss_coef() const {
    return _ss_coef;
}

double HamiltonianParamsAfFo::get_tau_z_coef() const {
    return _hamiltonian_params_fo.get_tau_z_coef();
}

double HamiltonianParamsAfFo::get_tau_minus_coef() const {
    return _hamiltonian_params_fo.get_tau_minus_coef();
}

double HamiltonianParamsAfFo::get_Pzz_coef() const {
    return _hamiltonian_params_fo.get_Pzz_coef();
}

double HamiltonianParamsAfFo::get_Pxz_coef() const {
    return _hamiltonian_params_fo.get_Pxz_coef();
}

double HamiltonianParamsAfFo::get_Pxx_coef() const {
    return _hamiltonian_params_fo.get_Pxx_coef();
}

double HamiltonianParamsAfFo::get_ss_Pzz_coef() const {
    return _hamiltonian_params_ss_fo.get_Pzz_coef();
}

double HamiltonianParamsAfFo::get_ss_Pxz_coef() const {
    return _hamiltonian_params_ss_fo.get_Pxz_coef();
}

double HamiltonianParamsAfFo::get_ss_Pxx_coef() const {
    return _hamiltonian_params_ss_fo.get_Pxx_coef();
}

monostar_hamiltonians::HamiltonianParamsFo HamiltonianParamsAfFo::get_hamiltonian_params_fo() const{
    return _hamiltonian_params_fo;
}

monostar_hamiltonians::HamiltonianParamsFo HamiltonianParamsAfFo::get_hamiltonian_params_ss_fo() const{
    return _hamiltonian_params_ss_fo;
}

double HamiltonianParamsAfFo::get_site_energy(double theta) const {
    return hamiltonian_params_af_fo_to_classic_energy_function(*this).get_value(theta);
}

double HamiltonianParamsAfFo::get_site_energy_derivative(double theta) const {
    return hamiltonian_params_af_fo_to_classic_energy_function(*this).get_derivative_value(theta);
}

double HamiltonianParamsAfFo::get_site_energy_derivative2(double theta) const {
    return hamiltonian_params_af_fo_to_classic_energy_function(*this).get_derivative2_value(theta);
}

double HamiltonianParamsAfFo::get_site_energy_derivative3(double theta) const {
    return hamiltonian_params_af_fo_to_classic_energy_function(*this).get_derivative3_value(theta);
}

double HamiltonianParamsAfFo::get_site_energy_derivative4(double theta) const {
    return hamiltonian_params_af_fo_to_classic_energy_function(*this).get_derivative4_value(theta);
}

std::set<double> HamiltonianParamsAfFo::get_theta_opt() const {
    return hamiltonian_params_af_fo_to_classic_energy_function(*this).get_minimum_argument();
}

std::set<double> HamiltonianParamsAfFo::get_theta_opt_numerical() const {
    return hamiltonian_params_af_fo_to_classic_energy_function(*this).get_minimum_argument_numerical();
}

utility::Result<std::set<double>, monostar_hamiltonians::NoKnownAnalyticalSolutionError> HamiltonianParamsAfFo::get_theta_opt_analytical() const {
    return hamiltonian_params_af_fo_to_classic_energy_function(*this).get_minimum_argument_analytical();
}

std::string HamiltonianParamsAfFo::string_repr_in_orbital_operators() const {
    std::ostringstream ss;
    ss << std::showpos;
    assert(_hamiltonian_params_ss_fo.get_tau_z_coef() == 0.0);
    assert(_hamiltonian_params_ss_fo.get_tau_minus_coef() == 0.0);
    if (_hamiltonian_params_fo.get_tau_z_coef() != 0) {
        ss << _hamiltonian_params_fo.get_tau_z_coef() << "·Σᵢ[τᶻᵢ]";
    }
    if (_hamiltonian_params_fo.get_tau_minus_coef() != 0) {
        ss << _hamiltonian_params_fo.get_tau_minus_coef() << "·Σᵢ[τ⁻ᵢ]";
    }
    if (_hamiltonian_params_fo.get_Pzz_coef() != 0) {
        ss << _hamiltonian_params_fo.get_Pzz_coef() << "·Σᵢⱼ[PᶻᵢPᶻⱼ]";
    }
    if (_hamiltonian_params_fo.get_Pxz_coef() != 0) {
        ss << _hamiltonian_params_fo.get_Pxz_coef() << "·Σᵢⱼ[PˣᵢPᶻⱼ + PᶻᵢPˣⱼ]";
    }
    if (_hamiltonian_params_fo.get_Pxx_coef() != 0) {
        ss << _hamiltonian_params_fo.get_Pxx_coef() << "·Σᵢⱼ[PˣᵢPˣⱼ]";
    }
    if (_hamiltonian_params_ss_fo.get_Pzz_coef() != 0) {
        ss << _hamiltonian_params_ss_fo.get_Pzz_coef() << "·Σᵢⱼ[SᵢSⱼPᶻᵢPᶻⱼ]";
    }
    if (_hamiltonian_params_ss_fo.get_Pxz_coef() != 0) {
        ss << _hamiltonian_params_ss_fo.get_Pxz_coef() << "·Σᵢⱼ[SᵢSⱼ[PˣᵢPᶻⱼ + PᶻᵢPˣⱼ]]";
    }
    if (_hamiltonian_params_ss_fo.get_Pxx_coef() != 0) {
        ss << _hamiltonian_params_ss_fo.get_Pxx_coef() << "·Σᵢⱼ[SᵢSⱼPˣᵢPˣⱼ]";
    }
    if (_hamiltonian_params_fo.get_tau_z_coef() == 0 &&
        _hamiltonian_params_fo.get_tau_minus_coef() == 0 &&
        _hamiltonian_params_fo.get_Pzz_coef() == 0 &&
        _hamiltonian_params_fo.get_Pxz_coef() == 0 &&
        _hamiltonian_params_fo.get_Pxx_coef() == 0 &&
        _hamiltonian_params_ss_fo.get_Pzz_coef() == 0 &&
        _hamiltonian_params_ss_fo.get_Pxz_coef() == 0 &&
        _hamiltonian_params_ss_fo.get_Pxx_coef() == 0) {
        ss << 0.0;
    }
    return ss.str();
}

std::string HamiltonianParamsAfFo::string_repr_in_trigonometric_functions() const {
    double SS_average = -4.0; //TODO: make SS_average an arg.
    std::ostringstream ss;
    ss << std::showpos;
    ss << _hamiltonian_params_fo.string_repr_in_trigonometric_functions();
    ss << SS_average << "·[" << _hamiltonian_params_ss_fo.string_repr_in_orbital_operators() << "]";
    return ss.str();
}

}  // end of namespace so_hamiltonians
