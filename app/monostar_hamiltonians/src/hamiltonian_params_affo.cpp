#include <monostar_hamiltonians/hamiltonian_params_affo.hpp>
#include <monostar_hamiltonians/hamiltonian_params_affo_to_classic_energy_function.hpp>

#include <monostar_hamiltonians/hamiltonian_params_fo_helpers.hpp>
#include <monostar_hamiltonians/hamiltonian_params_fo_site_matrices.hpp>
#include <monostar_hamiltonians/hamiltonian_params_fo_to_classic_energy_function.hpp>

#include <sstream>
#include <iomanip>

// #######################################################################
// ## HamiltonianParamsAffo                                             ##
// #######################################################################

namespace monostar_hamiltonians {

//HamiltonianParamsAffo::Builder HamiltonianParamsAffo::Builder::set_s_coef(double s_coef) {
//    _s_coef = s_coef;
//    return *this;
//}

HamiltonianParamsAffo::Builder HamiltonianParamsAffo::Builder::set_ss_coef(double ss_coef) {
    _ss_coef = ss_coef;
    return *this;
}

HamiltonianParamsAffo::Builder HamiltonianParamsAffo::Builder::set_tau_z_coef(double tau_z_coef) {
    _tau_z_coef = tau_z_coef;
    return *this;
}

HamiltonianParamsAffo::Builder HamiltonianParamsAffo::Builder::set_tau_minus_coef(double tau_minus_coef) {
    _tau_munis_coef = tau_minus_coef;
    return *this;
}

HamiltonianParamsAffo::Builder HamiltonianParamsAffo::Builder::set_Pzz_coef(double Pzz_coef) {
    _Pzz_coef = Pzz_coef;
    return *this;
}

HamiltonianParamsAffo::Builder HamiltonianParamsAffo::Builder::set_Pxz_coef(double Pxz_coef) {
    _Pxz_coef = Pxz_coef;
    return *this;
}

HamiltonianParamsAffo::Builder HamiltonianParamsAffo::Builder::set_Pxx_coef(double Pxx_coef) {
    _Pxx_coef = Pxx_coef;
    return *this;
}

HamiltonianParamsAffo::Builder HamiltonianParamsAffo::Builder::set_ss_Pzz_coef(double ss_Pzz_coef) {
    _ss_Pzz_coef = ss_Pzz_coef;
    return *this;
}

HamiltonianParamsAffo::Builder HamiltonianParamsAffo::Builder::set_ss_Pxz_coef(double ss_Pxz_coef) {
    _ss_Pxz_coef = ss_Pxz_coef;
    return *this;
}

HamiltonianParamsAffo::Builder HamiltonianParamsAffo::Builder::set_ss_Pxx_coef(double ss_Pxx_coef) {
    _ss_Pxx_coef = ss_Pxx_coef;
    return *this;
}

HamiltonianParamsAffo HamiltonianParamsAffo::Builder::build() const {
    return HamiltonianParamsAffo(
        /*_s_coef,*/ _ss_coef,
        _tau_z_coef, _tau_munis_coef, _Pzz_coef, _Pxz_coef, _Pxx_coef,
        _ss_Pzz_coef, _ss_Pxz_coef, _ss_Pxx_coef);
}

// #######################################################################

HamiltonianParamsAffo::HamiltonianParamsAffo(
    /*double s_coef,*/ double ss_coef,
    double tau_z_coef, double tau_minus_coef,
    double Pzz_coef, double Pxz_coef, double Pxx_coef,
    double ss_Pzz_coef, double ss_Pxz_coef, double ss_Pxx_coef)
    : /*_s_coef(s_coef),*/
      _ss_coef(ss_coef),
      _hamiltonian_params_fo(HamiltonianParamsFo::Builder()
                                 .set_tau_z_coef(tau_z_coef)
                                 .set_tau_minus_coef(tau_minus_coef)
                                 .set_Pzz_coef(Pzz_coef)
                                 .set_Pxz_coef(Pxz_coef)
                                 .set_Pxx_coef(Pxx_coef)
                                 .build()),
      _hamiltonian_params_ss_fo(HamiltonianParamsFo::Builder()
                                    .set_tau_z_coef(0.0)
                                    .set_tau_minus_coef(0.0)
                                    .set_Pzz_coef(ss_Pzz_coef)
                                    .set_Pxz_coef(ss_Pxz_coef)
                                    .set_Pxx_coef(ss_Pxx_coef)
                                    .build()) {
}

//double HamiltonianParamsAffo::get_s_coef() const {
//    return _s_coef;
//}

double HamiltonianParamsAffo::get_ss_coef() const {
    return _ss_coef;
}

double HamiltonianParamsAffo::get_tau_z_coef() const {
    return _hamiltonian_params_fo.get_tau_z_coef();
}

double HamiltonianParamsAffo::get_tau_minus_coef() const {
    return _hamiltonian_params_fo.get_tau_minus_coef();
}

double HamiltonianParamsAffo::get_Pzz_coef() const {
    return _hamiltonian_params_fo.get_Pzz_coef();
}

double HamiltonianParamsAffo::get_Pxz_coef() const {
    return _hamiltonian_params_fo.get_Pxz_coef();
}

double HamiltonianParamsAffo::get_Pxx_coef() const {
    return _hamiltonian_params_fo.get_Pxx_coef();
}

double HamiltonianParamsAffo::get_ss_Pzz_coef() const {
    return _hamiltonian_params_ss_fo.get_Pzz_coef();
}

double HamiltonianParamsAffo::get_ss_Pxz_coef() const {
    return _hamiltonian_params_ss_fo.get_Pxz_coef();
}

double HamiltonianParamsAffo::get_ss_Pxx_coef() const {
    return _hamiltonian_params_ss_fo.get_Pxx_coef();
}

HamiltonianParamsFo
HamiltonianParamsAffo::average_out_spins_1(double average_s) const {
    return average_out_spins_12(average_s, -average_s * average_s);
}

HamiltonianParamsFo
HamiltonianParamsAffo::average_out_spins_12([[maybe_unused]] double average_s, double average_ss) const {
    const double free =
        /*+get_s_coef() * average_s*/ +get_ss_coef() * average_ss;
    return HamiltonianParamsFo::Builder()
        .set_tau_z_coef(get_tau_z_coef())
        .set_tau_minus_coef(get_tau_minus_coef())
        .set_Pzz_coef(get_Pzz_coef() + average_ss * get_ss_Pzz_coef())
        .set_Pxz_coef(get_Pxz_coef() + average_ss * get_ss_Pxz_coef())
        .set_Pxx_coef(get_Pxx_coef() + average_ss * get_ss_Pxx_coef())
        .set_free_coef(free)
        .build();
}

HamiltonianParamsAfFm
HamiltonianParamsAffo::average_out_orbitals_1(double theta) const {
        const double average_tau_minus = std::real(OneSiteOrbitalMatrices::get_tau_minus_in_ge_basis(theta)(0, 0));
        const double average_tau_z = std::real(OneSiteOrbitalMatrices::get_tau_z_in_ge_basis(theta)(0, 0));
        const double average_Pzz = std::real(TwoSitesOrbitalMatrices::get_P_zz_in_ge_basis(theta)(0, 0));
        const double average_Pzx_sum_P_xz = std::real(TwoSitesOrbitalMatrices::get_P_zx_sum_P_xz_in_ge_basis(theta)(0, 0));
        const double average_Pxx = std::real(TwoSitesOrbitalMatrices::get_P_xx_in_ge_basis(theta)(0, 0));
        return average_out_orbitals_12(
            average_tau_minus, average_tau_z,
            average_Pzz, average_Pzx_sum_P_xz, average_Pxx);
}

HamiltonianParamsAfFm
HamiltonianParamsAffo::average_out_orbitals_12(
    double average_tau_minus, double average_tau_z,
    double average_Pzz, double average_Pzx_sum_P_xz, double average_Pxx) const {
        const double free =
            +get_tau_minus_coef() * average_tau_minus + get_tau_z_coef() * average_tau_z + get_Pzz_coef() * average_Pzz + get_Pxz_coef() * average_Pzx_sum_P_xz + get_Pxx_coef() * average_Pxx;
        const double B =
            /*+get_s_coef()*/ 0.0;
        const double J =
            +get_ss_coef() + get_ss_Pzz_coef() * average_Pzz + get_ss_Pxz_coef() * average_Pzx_sum_P_xz + get_ss_Pxx_coef() * average_Pxx;
        return HamiltonianParamsAfFm::Builder()
            .set_B(B)
            .set_J_classical(J)
            .set_J_quantum(J)
            .set_free(free)
            .build();
}

HamiltonianParamsFo HamiltonianParamsAffo::get_hamiltonian_params_fo() const {
    return _hamiltonian_params_fo;
}

HamiltonianParamsFo HamiltonianParamsAffo::get_hamiltonian_params_ss_fo() const {
    return _hamiltonian_params_ss_fo;
}

double HamiltonianParamsAffo::get_site_energy(double theta, double average_ss) const {
    return hamiltonian_params_affo_to_classic_energy_function_of_theta(*this, average_ss).get_value(theta);
}

double HamiltonianParamsAffo::get_site_energy_derivative(double theta, double average_ss) const {
    return hamiltonian_params_affo_to_classic_energy_function_of_theta(*this, average_ss).get_derivative_value(theta);
}

double HamiltonianParamsAffo::get_site_energy_derivative2(double theta, double average_ss) const {
    return hamiltonian_params_affo_to_classic_energy_function_of_theta(*this, average_ss).get_derivative2_value(theta);
}

double HamiltonianParamsAffo::get_site_energy_derivative3(double theta, double average_ss) const {
    return hamiltonian_params_affo_to_classic_energy_function_of_theta(*this, average_ss).get_derivative3_value(theta);
}

double HamiltonianParamsAffo::get_site_energy_derivative4(double theta, double average_ss) const {
    return hamiltonian_params_affo_to_classic_energy_function_of_theta(*this, average_ss).get_derivative4_value(theta);
}

std::set<double> HamiltonianParamsAffo::get_theta_opt(double average_ss) const {
    return hamiltonian_params_affo_to_classic_energy_function_of_theta(*this, average_ss).get_minimum_argument();
}

std::set<double> HamiltonianParamsAffo::get_theta_opt_numerical(double average_ss) const {
    return hamiltonian_params_affo_to_classic_energy_function_of_theta(*this, average_ss).get_minimum_argument_numerical();
}

utility::Result<std::set<double>, NoKnownAnalyticalSolutionError> HamiltonianParamsAffo::get_theta_opt_analytical(double average_ss) const {
    return hamiltonian_params_affo_to_classic_energy_function_of_theta(*this, average_ss).get_minimum_argument_analytical();
}

std::string HamiltonianParamsAffo::string_repr_in_orbital_operators() const {
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

std::string HamiltonianParamsAffo::string_repr_in_trigonometric_functions() const {
    double average_ss = -4.0;  //TODO: make average_ss an arg.
    std::ostringstream ss;
    ss << std::showpos;
    ss << _hamiltonian_params_fo.string_repr_in_trigonometric_functions();
    ss << average_ss << "·[" << _hamiltonian_params_ss_fo.string_repr_in_orbital_operators() << "]";
    return ss.str();
}

}  // end of namespace monostar_hamiltonians
