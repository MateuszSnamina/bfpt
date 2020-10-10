#include<monostar_hamiltonians/hamiltonian_params_fo.hpp>
#include<monostar_hamiltonians/hamiltonian_params_fo_helpers.hpp>

// #######################################################################
// ## hamiltonian_params_fo_to_classic_energy_function                  ##
// #######################################################################

namespace {

using namespace monostar_hamiltonians;

AcosPlusBsinPlusCsqcosPlusZ hamiltonian_params_fo_to_classic_energy_function(HamiltonianParamsFo params) {
    // E(θ) = + A*cos(θ) + B*sin(θ)  + C*cos⁴(θ/2) + 2*D*cos²(θ/2)sin²(θ/2) + E*sin⁴(θ/2)
    //      = + A*cos(θ) + B*sin(θ)  + C[½+½cos(θ)]² + 2*D[½*sin(θ)]² + E[½-½cos(θ)]²
    //      = + A*cos(θ) + B*sin(θ)  + ¼C[1+cos(θ)]² + ½D[sin(θ)]² + ¼E[1-cos(θ)]²
    //      = + A*cos(θ) + B*sin(θ)  + (¼C+¼E) + ½(C-E)*cos(θ) + ¼(C+E)*cos²(θ) + ½D - ½D*cos²(θ)
    //      =  [A+½(C-E)]*cos(θ) - B*sin(θ) + ¼(C+E-2D)*cos²(θ) + ¼(C+E+2D)
    // where: A ≡ +tau_z_coef,
    //        B ≡ -tau_minus_coef,
    //        C ≡ +Pzz_coef,
    //        D ≡ +Pxz_coef,
    //        E ≡ +Pxx_coef.
    const double cos_coef = params.get_tau_z_coef() + 0.5 * ( + params.get_Pzz_coef() - params.get_Pxx_coef());
    const double sin_coef = -params.get_tau_minus_coef();
    const double sqcos_coef = 0.25 * (params.get_Pzz_coef() + params.get_Pxx_coef() - 2 * params.get_Pxz_coef());
    const double free_coef = 0.25 * (params.get_Pzz_coef() + params.get_Pxx_coef() + 2 * params.get_Pxz_coef());
    return AcosPlusBsinPlusCsqcosPlusZ::Builder()
            .set_cos_coef(cos_coef)
            .set_sin_coef(sin_coef)
            .set_sqcos_coef(sqcos_coef)
            .set_free_coef(free_coef)
            .build();
}

} // end of namespace

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

HamiltonianParamsFo HamiltonianParamsFo::Builder::build() const {
    return HamiltonianParamsFo(_tau_z_coef, _tau_munis_coef, _Pzz_coef, _Pxz_coef, _Pxx_coef);
}

// #######################################################################

HamiltonianParamsFo::HamiltonianParamsFo(double tau_z_coef, double tau_minus_coef, double Pzz_coef, double Pxz_coef, double Pxx_coef) :
    _tau_z_coef(tau_z_coef),
    _tau_minus_coef(tau_minus_coef),
    _Pzz_coef(Pzz_coef),
    _Pxz_coef(Pxz_coef),
    _Pxx_coef(Pxx_coef) {
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

double HamiltonianParamsFo::get_site_energy(double theta) const {
    return hamiltonian_params_fo_to_classic_energy_function(*this).get_value(theta);
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

} // end of namespace monostar_hamiltonians
