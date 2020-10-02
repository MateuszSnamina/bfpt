#pragma once

#include<set>

// #######################################################################
// ## AcosPlusBsinPlusZ                                                 ##
// #######################################################################

/*
 * f(θ) = A*cos(θ) + B*sin(θ) + Z
 *      = H * [A/H * cos(θ) + B/H * sin(θ)]           // where: H ≡ √(A²+B²), valid if H≠0
 *      = H * [(-y) * cos(θ) + (x) * sin(θ)]          // where: y ≡ -A/H, x ≡ +B/H
 *      = H * [-sin(θ₀) * cos(θ) + cos(θ₀) * sin(θ)]  // where: θ₀ ≡ atan2(y, x);
 *      = H * sin(θ-θ₀)
 * f(θ) = A*cos(θ) + B*sin(θ) + Z
 *      = H * [A/H * cos(θ) + B/H * sin(θ)]           // where: H ≡ √(A²+B²), valid if H≠0
 *      = -H * [x * cos(θ) + y * sin(θ)]               // where: x ≡ -A/H, y ≡ -B/H
 *      = -H * [cos(θ₀) * cos(θ) + sin(θ₀) * sin(θ)]   // where: θ₀ ≡ atan2(y, x);
 *      = -H * cos(θ-θ₀)
 */

class AcosPlusBsinPlusZ {
public:
    class Builder {
    public:
        Builder set_cos_coef(double);
        Builder set_sin_coef(double);
        Builder set_free_coef(double);
        AcosPlusBsinPlusZ build() const;
    private:
        double _cos_coef = 0.0;
        double _sin_coef = 0.0;
        double _free_coef = 0.0;
    };
    double get_cos_coef() const;
    double get_sin_coef() const;
    double get_free_coef() const;
    double get_value(double phi) const ;
    std::set<double> get_minimum_argument() const;
private:
    AcosPlusBsinPlusZ(double cos_coef, double sin_coef, double free_coef);
    bool is_degenerated_to_const_function() const;
    const double _cos_coef;
    const double _sin_coef;
    const double _free_coef;
};

// #######################################################################
// ## BsinPlusCsqcosPlusZ                                               ##
// #######################################################################

/*
 * f(θ) = B*sin(θ) + C*cos²(θ) + Z
 * df/dθ = B*cos(θ) - 2*C*cos(θ)*sin(θ)
 *       = -2C * [-B'*cos(θ) + cos(θ)*sin(θ)]     where B'≡+B/2C if C≠0
 *       = -2C * [-B' + sin(θ)] * cos(θ)
 * d²f/dθ² = 2C * [-B' + sin(θ)] * sin(θ) - 2C * cos²(θ)
 * Case when B' < -1:
 * θ₁ = π/2
 * θ₂ = 3π/2
 * d²f/dθ²[θ₁] = +2C * [-B' + 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1-B') = +sgn(C)
 * d²f/dθ²[θ₂] = -2C * [-B' - 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1+B') = -sgn(C)
 * θ₁ is a minimum if C > 0.
 * Case when B' = -1:
 * θ₁ = π/2
 * θ₂ = -π/2
 * d²f/dθ²[θ₁] = +2C * [-B' + 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1-B') = +sgn(C)
 * d²f/dθ²[θ₂] = -2C * [-B' - 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1+B') = 0
 * ⇒ θ₁ is a minimum if C > 0.
 * Case when -1 < B' < +1:
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
 * Case when B' = +1:
 * θ₁ = π/2
 * θ₂ = 3π/2
 * d²f/dθ²[θ₁] = +2C * [-B' + 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1-B') = 0
 * d²f/dθ²[θ₂] = -2C * [-B' - 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1+B') = +sgn(C)
 * ⇒ θ₂ is minimum if C > 0.
 * Case when B' > 1:
 * θ₁ = π/2
 * θ₂ = -π/2
 * d²f/dθ²[θ₁] = +2C * [-B' + 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1-B') = -sgn(C)
 * d²f/dθ²[θ₂] = -2C * [-B' - 1] ⇒ sgn(d²f/dθ²[θ₁]) = +sgn(C)sgn(1+B') = +sgn(C)
 * ⇒ θ₂ is minimum if C > 0.
 */

class BsinPlusCsqcosPlusZ {
public:
    class Builder {
    public:
        Builder set_sin_coef(double);
        Builder set_sqcos_coef(double);
        Builder set_free_coef(double);
        BsinPlusCsqcosPlusZ build() const;
    private:
        double _sin_coef = 0.0;
        double _sqcos_coef = 0.0;
        double _free_coef = 0.0;
    };
    double get_sin_coef() const;
    double get_sqcos_coef() const;
    double get_free_coef() const;
    double get_value(double phi) const;
    std::set<double> get_minimum_argument() const;
private:
    BsinPlusCsqcosPlusZ(double sin_coef, double sqcos_coef, double free_coef);
    const double _sin_coef;
    const double _sqcos_coef;
    const double _free_coef;
    bool is_degenerated_to_const_function() const;
    bool is_degenerated_to_sqcos_coef_equal_to_zero_case() const;
    std::set<double> get_minimum_argument_analytical_when_sqcos_coef_is_zero() const;
    std::set<double> get_minimum_argument_analytical_when_sqcos_coef_is_not_zero() const;
};

//// #######################################################################
//// ## AcosPlusCsqcosPlusZ                                               ##
//// #######################################################################

///*
// * f(θ) = A*cos(θ) + C*cos²(θ) + Z
// * df/dθ = -A*sin(θ) - 2*C*cos(θ)*sin(θ)
// *       = -2C * [-A'*sin(θ) +cos(θ)*sin(θ)]     where A'≡-A/2C if C≠0
// *       = -2C * [sin(θ)] * [-A' + cos(θ)] + 2*C*A'*B'
// */

//class AcosPlusCsqcosPlusZ {
//public:
//    class Builder {
//    public:
//        Builder set_cos_coef(double);
//        Builder set_sqcos_coef(double);
//        Builder set_free_coef(double);
//        AcosPlusCsqcosPlusZ build() const;
//    private:
//        double _cos_coef = 0.0;
//        double _sqcos_coef = 0.0;
//        double _free_coef = 0.0;
//    };
//    double get_cos_coef() const;
//    double get_sqcos_coef() const;
//    double get_free_coef() const;
//    double get_value(double phi) const;
//    std::set<double> get_minimum_argument() const;
//private:
//    AcosPlusCsqcosPlusZ(double cos_coef, double sqcos_coef, double free_coef);
//    const double _cos_coef;
//    const double _sqcos_coef;
//    const double _free_coef;
//};
