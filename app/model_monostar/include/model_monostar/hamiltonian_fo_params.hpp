#pragma once

#include<optional>

// #######################################################################
// ## HamiltonianFoParams                                               ##
// #######################################################################

class HamiltonianFoParams {
public:
    double get_tau_minus_coef() const;
    double get_tau_z_coef() const;
    double get_Pzz_coef() const;
    double get_Pxz_coef() const;
    double get_Pxx_coef() const;
    class Builder {
    public:
        Builder set_tau_minus_coef(double);
        Builder set_tau_z_coef(double);
        Builder set_Pzz_coef(double);
        Builder set_Pxz_coef(double);
        Builder set_Pxx_coef(double);
        HamiltonianFoParams build() const;
    private:
        double _tau_munis_coef = 0.0;
        double _tau_z_coef = 0.0;
        double _Pzz_coef = 0.0;
        double _Pxz_coef = 0.0;
        double _Pxx_coef = 0.0;
    };
    friend HamiltonianFoParams Builder::build() const;
    HamiltonianFoParams() = default;
    double get_site_energy(double theta);
    double get_theta_opt();
private:
    HamiltonianFoParams(double tau_minus_coef, double tau_z_coef, double Pzz_coef, double Pxz_coef, double Pxx_coef);
    double _tau_minus_coef;
    double _tau_z_coef;
    double _Pzz_coef;
    double _Pxz_coef;
    double _Pxx_coef;
};
