#include<model_monostar/hamiltonian_fo_params.hpp>
#include<model_monostar/hamiltonian_fo_params_helpers.hpp>

// #######################################################################
// ## hamiltonian_fo_params_to_classic_energy_function                  ##
// #######################################################################

namespace {

AcosPlucBsinPlusCsqcosPlusZ hamiltonian_fo_params_to_classic_energy_function(HamiltonianFoParams params) {
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
    return AcosPlucBsinPlusCsqcosPlusZ::Builder()
            .set_cos_coef(cos_coef)
            .set_sin_coef(sin_coef)
            .set_sqcos_coef(sqcos_coef)
            .set_free_coef(free_coef)
            .build();
}

} // end of namespace

// #######################################################################
// ## HamiltonianFoParams                                               ##
// #######################################################################

HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_tau_z_coef(double tau_z_coef) {
    _tau_z_coef = tau_z_coef;
    return *this;
}

HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_tau_minus_coef(double tau_minus_coef) {
    _tau_munis_coef = tau_minus_coef;
    return *this;
}

HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_Pzz_coef(double Pzz_coef) {
    _Pzz_coef = Pzz_coef;
    return *this;
}

HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_Pxz_coef(double Pxz_coef) {
    _Pxz_coef = Pxz_coef;
    return *this;
}

HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_Pxx_coef(double Pxx_coef) {
    _Pxx_coef = Pxx_coef;
    return *this;
}

HamiltonianFoParams HamiltonianFoParams::Builder::build() const {
    return HamiltonianFoParams(_tau_z_coef, _tau_munis_coef, _Pzz_coef, _Pxz_coef, _Pxx_coef);
}

// #######################################################################

HamiltonianFoParams::HamiltonianFoParams(double tau_z_coef, double tau_minus_coef, double Pzz_coef, double Pxz_coef, double Pxx_coef) :
    _tau_z_coef(tau_z_coef),
    _tau_minus_coef(tau_minus_coef),
    _Pzz_coef(Pzz_coef),
    _Pxz_coef(Pxz_coef),
    _Pxx_coef(Pxx_coef) {
}

double HamiltonianFoParams::get_tau_z_coef() const {
    return _tau_z_coef;
}

double HamiltonianFoParams::get_tau_minus_coef() const {
    return _tau_minus_coef;
}

double HamiltonianFoParams::get_Pzz_coef() const {
    return _Pzz_coef;
}

double HamiltonianFoParams::get_Pxz_coef() const {
    return _Pxz_coef;
}

double HamiltonianFoParams::get_Pxx_coef() const {
    return _Pxx_coef;
}

double HamiltonianFoParams::get_site_energy(double theta) const {
    return hamiltonian_fo_params_to_classic_energy_function(*this).get_value(theta);
}

std::set<double> HamiltonianFoParams::get_theta_opt() const {
    return hamiltonian_fo_params_to_classic_energy_function(*this).get_minimum_argument();
}
