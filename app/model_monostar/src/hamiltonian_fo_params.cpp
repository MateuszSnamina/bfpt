#include<model_monostar/hamiltonian_fo_params.hpp>

#include<cmath>
#include<limits>
#include<cassert>

// https://stackoverflow.com/questions/1727881/how-to-use-the-pi-constant-in-c
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// #######################################################################
// ## almost_equal                                                      ##
// #######################################################################

// Ref: https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon

namespace {

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp = 100) {
    return std::fabs(x-y) <= std::numeric_limits<T>::epsilon() * std::fabs(x+y) * ulp
            || std::fabs(x-y) < std::numeric_limits<T>::min();
}

}

// #######################################################################
// ## AcosPlucBsinPlusZ                                                 ##
// #######################################################################

/*
 * f(θ) = A*cos(θ) + B*sin(θ) + Z
 *      = H * [A/H * cos(θ) + B/H * sin(θ)]           // where: H ≡ √(A²+B²), valid if H≠0
 *      = H * [(-y) * cos(θ) + (x) * sin(θ)]          // where: y ≡ -A/H, x ≡ +B/H
 *      = H * [-sin(θ₀) * cos(θ) + cos(θ₀) * sin(θ)]  // where: θ₀ ≡ atan2(y, x);
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
    bool is_degenerated_to_const_function() const;
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
    if (is_degenerated_to_const_function()) {
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

bool AcosPlucBsinPlusZ::is_degenerated_to_const_function() const {
    const double h = std::hypot(_cos_coef, _sin_coef);
    return h < 100 * std::numeric_limits<double>::epsilon();
}

} // end of namespace

// #######################################################################
// ## AcosPlucBsinPlusCsqcosPlusZ                                       ##
// #######################################################################

/*
 * f(θ) = A*cos(θ) + B*sin(θ) + C*cos²(θ) + Z
 * df/dθ = -A*sin(θ) + B*cos(θ) - 2*C*cos(θ)*sin(θ)
 *       = -2C * [-A'*sin(θ) -B'*cos(θ) +cos(θ)*sin(θ)]     where A'≡-A/2C, B'≡+B/2C if C≠0
 *       = -2C * [-B' + sin(θ)] * [-A' + cos(θ)] + 2*C*A'*B'
 *
 * ******** Case when A=0 ********
 * df/dθ = -2C * [-B' + sin(θ)] * cos(θ)
 * d²f/dθ² = 2C * [-B' + sin(θ)] * sin(θ) - 2C * cos²(θ)
 * Case when A=0 ⋀ B' < -1:
 * θ₁ = π/2
 * θ₂ = 3π/2
 * d²f/dθ²[θ₁] = +2C * [-B' + 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1-B') = +sgn(C)
 * d²f/dθ²[θ₂] = -2C * [-B' - 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1+B') = -sgn(C)
 * θ₁ is a minimum if C > 0.
 * Case when A=0 ⋀ B' = -1:
 * θ₁ = π/2
 * θ₂ = 3π/2
 * d²f/dθ²[θ₁] = +2C * [-B' + 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1-B') = +sgn(C)
 * d²f/dθ²[θ₂] = -2C * [-B' - 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1+B') = 0
 * ⇒ θ₁ is a minimum if C > 0.
 * Case when A=0 ⋀ -1 < B' < +1:
 * θ₁ = π/2
 * θ₂ = 3π/2
 * θ₃ = arcsin(B')    // θ₃ is in -π/2..+π/2
 * θ₄ = π-arcsin(B')  // as: sin(π-θ)=sin(θ)
 * d²f/dθ²[θ₁] = +2C * [-B' + 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1-B') = +sgn(C)
 * d²f/dθ²[θ₂] = -2C * [-B' - 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1+B') = +sgn(C)
 * d²f/dθ²[θ₃] = -2C * cos²(θ₃) ⇒ sgn(d²f/dθ²[θ₃]) =  -sgn(C)
 * d²f/dθ²[θ₄] = -2C * cos²(θ₄) ⇒ sgn(d²f/dθ²[θ₄]) =  -sgn(C)
 * Note that: f(θ₁)= +B + Z
 * Note that: f(θ₂)= -B + Z
 * Note that: f(θ₃)=f(θ₄) as sin(θ₃)=sin(π-θ₃)=sin(θ₄) and cos(θ₃)=-cos(π-θ₃)=-cos(θ₄)
 * ⇒ θ₁ and θ₂ are local minimum if C > 0.
 * ⇒ θ₁ is global  minimum if C > 0 ⋀ B < 0.
 * ⇒ θ₂ is global  minimum if C > 0 ⋀ B > 0.
 * ⇒ both θ₁ and θ₂ are global minimum if C > 0 ⋀ B = 0.
 * ⇒ both θ₃ and θ₄ are global minimum if C < 0.
 * Case when A=0 ⋀ B' = +1:
 * θ₁ = π/2
 * θ₂ = 3π/2
 * d²f/dθ²[θ₁] = +2C * [-B' + 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1-B') = 0
 * d²f/dθ²[θ₂] = -2C * [-B' - 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1+B') = +sgn(C)
 * ⇒ θ₂ is minimum if C > 0.
 * Case when A=0 ⋀ B' > 1:
 * θ₁ = π/2
 * θ₂ = 3π/2
 * d²f/dθ²[θ₁] = +2C * [-B' + 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1-B') = -sgn(C)
 * d²f/dθ²[θ₂] = -2C * [-B' - 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1+B') = +sgn(C)
 * ⇒ θ₂ is minimum if C > 0.
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
    bool is_degenerated_to_const_function() const;
    bool is_degenerated_to_cos_coef_equal_to_zero_case() const;
    bool is_degenerated_to_sin_coef_equal_to_zero_case() const;
    bool is_degenerated_to_sqcos_coef_equal_to_zero_case() const;
    std::optional<double> get_minimum_argument_when_cos_coef_is_zero() const;
    std::optional<double> get_minimum_argument_when_sin_coef_is_zero() const;
    std::optional<double> get_minimum_argument_when_sqcos_coef_is_zero() const;
    std::optional<double> get_minimum_not_analitycal() const;
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
    if (is_degenerated_to_const_function()) {
        return std::nullopt;
    }
    if (is_degenerated_to_cos_coef_equal_to_zero_case()) {
        return get_minimum_argument_when_cos_coef_is_zero();
    }
    if (is_degenerated_to_sin_coef_equal_to_zero_case()) {
        return get_minimum_argument_when_sin_coef_is_zero();
    }
    if (is_degenerated_to_sqcos_coef_equal_to_zero_case()) {
        return get_minimum_argument_when_sqcos_coef_is_zero();
    }
    return get_minimum_not_analitycal();
}

std::optional<double> AcosPlucBsinPlusCsqcosPlusZ::get_minimum_argument_when_cos_coef_is_zero() const {
    assert(!is_degenerated_to_const_function());
    assert(almost_equal(_cos_coef, 0.0));
    assert(!almost_equal(_sqcos_coef, 0.0));
    if (_sin_coef <= -std::abs(2 * _sqcos_coef)) {
        return +M_PI/2;
    } else if (_sin_coef >= +std::abs(2 * _sqcos_coef)) {
        return -M_PI/2;
    } else {
        assert(std::abs(_sin_coef) < std::abs(2 * _sqcos_coef));
        if (_sqcos_coef > 0) {
            [[maybe_unused]] const double theta_1 = +M_PI/2;
            [[maybe_unused]] const double theta_2 = -M_PI/2;
            [[maybe_unused]] const double f_theta_1 = get_value(theta_1);
            [[maybe_unused]] const double f_theta_2 = get_value(theta_2);
            assert(f_theta_1 == +_sin_coef + _free_coef);
            assert(f_theta_2 == -_sin_coef + _free_coef);
            if (_sin_coef < 0) {
                assert(f_theta_1 <= f_theta_2 || almost_equal(f_theta_1, f_theta_2));
                return +M_PI/2;
            } else if (_sin_coef > 0) {
                assert(f_theta_2 <= f_theta_1 || almost_equal(f_theta_1, f_theta_2));
                return -M_PI/2;
            } else {
                assert(_sin_coef == 0);
                assert(almost_equal(f_theta_1, f_theta_2));
                return +M_PI/2;
                //return {-M_PI/2, +M_PI/2}; //TODO retrun both global minima.
            }
        } else {
            assert(std::abs(_sin_coef) < std::abs(2 * _sqcos_coef));
            assert(_sqcos_coef < 0);
            const double B_prim = _sin_coef/ (2 * _sqcos_coef);
            const double theta_3 = std::asin(B_prim);
            [[maybe_unused]] const double theta_4 = M_PI - std::asin(B_prim);
            [[maybe_unused]] const double f_theta_3 = get_value(theta_3);
            [[maybe_unused]] const double f_theta_4 = get_value(theta_4);
            assert(almost_equal(f_theta_3, f_theta_4));
            return f_theta_3;
            //return {f_theta_3, f_theta_3}; //TODO retrun both global minima.
        }
    }
}


std::optional<double> AcosPlucBsinPlusCsqcosPlusZ::get_minimum_argument_when_sin_coef_is_zero() const {
    assert(!is_degenerated_to_const_function());
    assert(almost_equal(_sin_coef, 0.0));
    assert(!almost_equal(_sqcos_coef, 0.0));
    assert(false);
    //TODO implement
    return 0.0/0.0;
}

std::optional<double> AcosPlucBsinPlusCsqcosPlusZ::get_minimum_argument_when_sqcos_coef_is_zero() const {
    return AcosPlucBsinPlusZ::Builder()
            .set_cos_coef(_cos_coef)
            .set_sin_coef(_sin_coef)
            .set_free_coef(_free_coef)
            .build()
            .get_minimum_argument();
}

std::optional<double> AcosPlucBsinPlusCsqcosPlusZ::get_minimum_not_analitycal() const {
    assert(false);
    //TODO implement
    return 0.0/0.0;
}

bool AcosPlucBsinPlusCsqcosPlusZ::is_degenerated_to_const_function() const {
    const double h = std::sqrt(_cos_coef * _cos_coef + _sqcos_coef * _sqcos_coef + _sin_coef * _sin_coef);
    return h < 100 * std::numeric_limits<double>::epsilon();
}


bool AcosPlucBsinPlusCsqcosPlusZ::is_degenerated_to_cos_coef_equal_to_zero_case() const {
    const double N = std::hypot(_sqcos_coef, _sin_coef);
    return std::abs(_cos_coef) < 100 * N * std::numeric_limits<double>::epsilon();
}

bool AcosPlucBsinPlusCsqcosPlusZ::is_degenerated_to_sin_coef_equal_to_zero_case() const {
    const double N = std::hypot(_cos_coef, _sqcos_coef);
    return std::abs(_sin_coef) < 100 * N * std::numeric_limits<double>::epsilon();
}

bool AcosPlucBsinPlusCsqcosPlusZ::is_degenerated_to_sqcos_coef_equal_to_zero_case() const {
    const double N = std::hypot(_cos_coef, _sin_coef);
    return std::abs(_sqcos_coef) < 100 * N * std::numeric_limits<double>::epsilon();
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
