#include<model_monostar/hamiltonian_fo_params.hpp>

#include<cmath>
#include<limits>
#include<cassert>

// https://stackoverflow.com/questions/1727881/how-to-use-the-pi-constant-in-c
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// #######################################################################
// ## AcosPlucBsinPlusZ                                                 ##
// #######################################################################

/*
 * f(θ) = A*cos(θ) + B*sin(θ) + Z
 *      = H * [A/H * cos(θ) + B/H * sin(θ)]           // where: H ≡ √(A²+B²), valid if H != 0
 *      = H * [(-y) * cos(θ) + (x) * sin(θ)]          // where: y = -A/H, x = +B/H
 *      = H * [-sin(θ₀) * cos(θ) + cos(θ₀) * sin(θ)]  // where: θ₀ = atan2(y, x);
 *      = H * sin(θ-θ₀)
 */

namespace {

class AcosPlucBsinPlusZ {
public:
    class Builder {
    public:
        Builder set_cos_coef(double);
        Builder set_sin_coef(double);
        Builder set_free_coef(double);
        AcosPlucBsinPlusZ build() const;
    private:
        double _cos_coef = 0.0;
        double _sin_coef = 0.0;
        double _free_coef = 0.0;
    };
    double get_cos_coef() const;
    double get_sin_coef() const;
    double get_free_coef() const;
    double get_value(double phi) const ;
    std::optional<double> get_minimum_argument() const;
private:
    AcosPlucBsinPlusZ(double cos_coef, double sin_coef, double free_coef);
    const double _cos_coef;
    const double _sin_coef;
    const double _free_coef;
};

AcosPlucBsinPlusZ::Builder AcosPlucBsinPlusZ::Builder::set_cos_coef(double cos_coef) {
    _cos_coef = cos_coef;
    return *this;
}

AcosPlucBsinPlusZ::Builder AcosPlucBsinPlusZ::Builder::set_sin_coef(double sin_coef) {
    _sin_coef = sin_coef;
    return *this;
}

AcosPlucBsinPlusZ::Builder AcosPlucBsinPlusZ::Builder::set_free_coef(double free_coef) {
    _free_coef = free_coef;
    return *this;
}

AcosPlucBsinPlusZ AcosPlucBsinPlusZ::Builder::build() const {
    return AcosPlucBsinPlusZ(_cos_coef, _sin_coef, _free_coef);
}

AcosPlucBsinPlusZ::AcosPlucBsinPlusZ(double cos_coef, double sin_coef, double free_coef) :
    _cos_coef(cos_coef),
    _sin_coef(sin_coef),
    _free_coef(free_coef) {
}

double AcosPlucBsinPlusZ::get_cos_coef() const {
    return _cos_coef;
}

double AcosPlucBsinPlusZ::get_sin_coef() const {
    return _sin_coef;
}

double AcosPlucBsinPlusZ::get_free_coef() const {
    return _free_coef;
}

double AcosPlucBsinPlusZ::get_value(double phi) const {
    return + _cos_coef * std::cos(phi)
            + _sin_coef * std::sin(phi)
            + _free_coef;
}

std::optional<double> AcosPlucBsinPlusZ::get_minimum_argument() const {
    const double h = std::hypot(_cos_coef, _sin_coef);
    if (h < std::numeric_limits<double>::epsilon()) {
        return std::nullopt;
    }
    const double x = + _sin_coef / h;
    const double y = - _cos_coef / h;
    const double phi_0 = std::atan2(y, x);
    if (std::abs(phi_0 + M_PI) < std::numeric_limits<double>::epsilon()) {
        return M_PI;
    } else {
        return phi_0;
    }
}

} // end of namespace

// #######################################################################
// ## AcosPlucBsinPlusCsqcosPlusZ                                       ##
// #######################################################################

/*
 * f(θ) = A*cos(θ) + B*sin(θ) + C*cos²(θ) + Z
 */

namespace {

class AcosPlucBsinPlusCsqcosPlusZ {
public:
    class Builder {
    public:
        Builder set_cos_coef(double);
        Builder set_sin_coef(double);
        Builder set_sqcos_coef(double);
        Builder set_free_coef(double);
        AcosPlucBsinPlusCsqcosPlusZ build() const;
    private:
        double _cos_coef = 0.0;
        double _sin_coef = 0.0;
        double _sqcos_coef = 0.0;
        double _free_coef = 0.0;
    };
    double get_cos_coef() const;
    double get_sin_coef() const;
    double get_sqcos_coef() const;
    double get_free_coef() const;
    double get_value(double phi) const ;
    std::optional<double> get_minimum_argument() const;
private:
    AcosPlucBsinPlusCsqcosPlusZ(double cos_coef, double sin_coef, double sqcos_coef, double free_coef);
    const double _cos_coef;
    const double _sin_coef;
    const double _sqcos_coef;
    const double _free_coef;
};

AcosPlucBsinPlusCsqcosPlusZ::Builder AcosPlucBsinPlusCsqcosPlusZ::Builder::set_cos_coef(double cos_coef) {
    _cos_coef = cos_coef;
    return *this;
}

AcosPlucBsinPlusCsqcosPlusZ::Builder AcosPlucBsinPlusCsqcosPlusZ::Builder::set_sin_coef(double sin_coef) {
    _sin_coef = sin_coef;
    return *this;
}

AcosPlucBsinPlusCsqcosPlusZ::Builder AcosPlucBsinPlusCsqcosPlusZ::Builder::set_sqcos_coef(double sqcos_coef) {
    _sqcos_coef = sqcos_coef;
    return *this;
}

AcosPlucBsinPlusCsqcosPlusZ::Builder AcosPlucBsinPlusCsqcosPlusZ::Builder::set_free_coef(double free_coef) {
    _free_coef = free_coef;
    return *this;
}

AcosPlucBsinPlusCsqcosPlusZ AcosPlucBsinPlusCsqcosPlusZ::Builder::build() const {
    return AcosPlucBsinPlusCsqcosPlusZ(_cos_coef, _sin_coef, _sqcos_coef, _free_coef);
}

AcosPlucBsinPlusCsqcosPlusZ::AcosPlucBsinPlusCsqcosPlusZ(double cos_coef, double sin_coef, double sqcos_coef, double free_coef) :
    _cos_coef(cos_coef),
    _sin_coef(sin_coef),
    _sqcos_coef(sqcos_coef),
    _free_coef(free_coef) {
}

double AcosPlucBsinPlusCsqcosPlusZ::get_cos_coef() const {
    return _cos_coef;
}

double AcosPlucBsinPlusCsqcosPlusZ::get_sin_coef() const {
    return _sin_coef;
}

double AcosPlucBsinPlusCsqcosPlusZ::get_sqcos_coef() const {
    return _sqcos_coef;
}

double AcosPlucBsinPlusCsqcosPlusZ::get_free_coef() const {
    return _free_coef;
}

double AcosPlucBsinPlusCsqcosPlusZ::get_value(double phi) const {
    return + _cos_coef * std::cos(phi)
            + _sin_coef * std::sin(phi)
            + _sqcos_coef * std::cos(phi) * std::cos(phi)
            + _free_coef;
}

std::optional<double> AcosPlucBsinPlusCsqcosPlusZ::get_minimum_argument() const {
    const double N = std::hypot(_cos_coef, _sin_coef);
    if (std::abs(_sqcos_coef) < 100 * N * std::numeric_limits<double>::epsilon()) {
        return AcosPlucBsinPlusZ::Builder()
                .set_cos_coef(_cos_coef)
                .set_sin_coef(_sin_coef)
                .set_free_coef(_free_coef)
                .build()
                .get_minimum_argument();
    } else {
        // TODO
        assert(false);
        return 0;
    }
}

} // end of namespace

// #######################################################################
// ## hamiltonian_fo_params_to_classic_energy_function                  ##
// #######################################################################

namespace {

AcosPlucBsinPlusCsqcosPlusZ hamiltonian_fo_params_to_classic_energy_function(HamiltonianFoParams params) {
    // E(θ) = + A*cos(θ) - B*sin(θ)  + C*cos⁴(θ/2) + 2*D*cos²(θ/2)sin²(θ/2) + E*sin⁴(θ/2)
    //      = + A*cos(θ) - B*sin(θ)  + C[½+½cos(θ)]² + 2*D[½*sin(θ)]² + E[½-½cos(θ)]²
    //      = + A*cos(θ) - B*sin(θ)  + ¼C[1+cos(θ)]² + ½D[sin(θ)]² + ¼E[1-cos(θ)]²
    //      = + A*cos(θ) - B*sin(θ)  + (¼C+¼E) + ½(C-E)*cos(θ) + ¼(C+E)*cos²(θ) + ½D - ½D*cos²(θ)
    //      =  [A+½(C-E)]*cos(θ) - B*sin(θ) + ¼(C+E-2D)*cos²(θ) + ¼(C+E+2D)
    // where: A ≡ tau_z_coef,
    //        B ≡ tau_minus_coef,
    //        C ≡ Pzz_coef,
    //        D ≡ Pxz_coef,
    //        E ≡ Pxx_coef.
    const double cos_coef = params.get_tau_z_coef() + 0.5 * ( + params.get_Pzz_coef() - params.get_Pxx_coef());
    const double sin_coef = params.get_tau_minus_coef();
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

std::optional<double> HamiltonianFoParams::get_theta_opt() const {
    return hamiltonian_fo_params_to_classic_energy_function(*this).get_minimum_argument();
}
