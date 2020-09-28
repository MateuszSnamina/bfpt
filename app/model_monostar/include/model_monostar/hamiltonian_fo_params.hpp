#pragma once

#include<cmath>

class HamiltonianFoParams {
public:
    double get_Pdelta_coef() const;
    double get_Pzz_coef() const;
    double get_Pxz_coef() const;
    double get_Pxx_coef() const;
    class Builder {
    public:
        Builder set_Pdelta_coef(double);
        Builder set_Pzz_coef(double);
        Builder set_Pxz_coef(double);
        Builder set_Pxx_coef(double);
        HamiltonianFoParams build() const;
    private:
        double _Pdelta_coef = 1.0;
        double _Pzz_coef = 1.0;
        double _Pxz_coef = 0.0;
        double _Pxx_coef = 0.0;
    };
    friend HamiltonianFoParams Builder::build() const;
    HamiltonianFoParams() = default;
    double get_site_energy(double theta);
    double get_theta_opt();
private:
    HamiltonianFoParams(double Bdelta, double Pzz_coef, double Pxz_coef, double Pxx_coef);
    double _Pdelta_coef = 1.0;
    double _Pzz_coef = 1.0;
    double _Pxz_coef = 0.0;
    double _Pxx_coef = 0.0;
};

inline
HamiltonianFoParams::HamiltonianFoParams(double Pdelta_coef, double Pzz_coef, double Pxz_coef, double Pxx_coef) :
    _Pdelta_coef(Pdelta_coef),
    _Pzz_coef(Pzz_coef),
    _Pxz_coef(Pxz_coef),
    _Pxx_coef(Pxx_coef) {
}

inline
double HamiltonianFoParams::get_Pdelta_coef() const{
    return _Pdelta_coef;
}

inline
double HamiltonianFoParams::get_Pzz_coef() const{
    return _Pzz_coef;
}

inline
double HamiltonianFoParams::get_Pxz_coef() const{
    return _Pxz_coef;
}

inline
double HamiltonianFoParams::get_Pxx_coef() const{
    return _Pxx_coef;
}

inline
double HamiltonianFoParams::get_site_energy(double theta) {
    const double s2 = std::sin(theta / 2);
    const double c2 = std::cos(theta / 2);
    const double s1 = std::sin(theta);
    const double Pdelta_gg = -s1;
    const double Px_gg = +s2 * s2;
    const double Pz_gg = +c2 * c2;
    const double H_1 = _Pdelta_coef * (Pdelta_gg);
    const double H_12zz = _Pzz_coef * (Pz_gg * Pz_gg);
    const double H_12xz = _Pxz_coef * (Px_gg * Pz_gg + Pz_gg * Px_gg);
    const double H_12xx = _Pxx_coef * (Px_gg * Px_gg);
    const double H_12 = H_12zz + H_12xz + H_12xx;
    return H_1 + H_12;
}

/*
 *
 * Function to optimize over θ:
 * E(θ) = - Asin(θ) + B*cos⁴(θ/2) + C*cos²(θ/2)sin²(θ/2) + D*sin⁴(θ/2)
 *      = - Asin(θ) + B[½+½cos(θ)]² + C[½*sin(θ)]² + D[½-½cos(θ)]²
 *      = - Asin(θ) + ¼B[1+cos(θ)]² + ¼C[sin(θ)]² + ¼D[1-cos(θ)]²
 *
 * dE/dθ = - Acos(θ) - ½B[1+cos(θ)]sin(θ) + ½C[sin(θ)]cos(θ) + ½D[1-cos(θ)]sin(θ)
 *       = - Acos(θ) + ½(-B+D)sin(θ) + ½(-B+C-D)sin(θ)cos(θ)
 * d²E/dθ² =
 *
 */

inline
double HamiltonianFoParams::get_theta_opt(){
    //assert(false);
    //TODO implement
}

inline
HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_Pdelta_coef(double Pdelta_coef) {
    _Pdelta_coef = Pdelta_coef;
    return *this;
}

inline
HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_Pzz_coef(double Pzz_coef) {
    _Pzz_coef = Pzz_coef;
    return *this;
}

inline
HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_Pxz_coef(double Pxz_coef) {
    _Pxz_coef = Pxz_coef;
    return *this;
}

inline
HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_Pxx_coef(double Pxx_coef) {
    _Pxx_coef = Pxx_coef;
    return *this;
}

inline
HamiltonianFoParams HamiltonianFoParams::Builder::build() const {
    return HamiltonianFoParams(_Pdelta_coef, _Pzz_coef, _Pxz_coef, _Pxx_coef);
}
