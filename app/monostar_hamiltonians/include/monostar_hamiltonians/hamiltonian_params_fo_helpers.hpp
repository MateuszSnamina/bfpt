#pragma once

#include <set>

#include <utility/result.hpp>

// #######################################################################
// ## NoKnownAnalyticalSolutionError                                    ##
// #######################################################################

namespace monostar_hamiltonians {

class NoKnownAnalyticalSolutionError : public std::domain_error {
   public:
    NoKnownAnalyticalSolutionError();
};

}  // end of namespace monostar_hamiltonians

// #######################################################################
// ## AcosPlusBsinPlusCsqcosPlusZ                                       ##
// #######################################################################

/*
 * f(θ) = A*cos(θ) + B*sin(θ) + C*cos²(θ) + Z
 * df/dθ = -A*sin(θ) + B*cos(θ) - 2*C*cos(θ)*sin(θ)
 *       = -2C * [-A'*sin(θ) -B'*cos(θ) +cos(θ)*sin(θ)]     where A'≡-A/2C, B'≡+B/2C if C≠0
 *       = -2C * [-B' + sin(θ)] * [-A' + cos(θ)] + 2*C*A'*B'
 */

namespace monostar_hamiltonians {

class AcosPlusBsinPlusCsqcosPlusZ {
   public:
    class Builder {
       public:
        Builder set_cos_coef(double);
        Builder set_sin_coef(double);
        Builder set_sqcos_coef(double);
        Builder set_free_coef(double);
        AcosPlusBsinPlusCsqcosPlusZ build() const;

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
    double get_value(double phi) const;
    double get_derivative_value(double phi) const;
    double get_derivative2_value(double phi) const;
    double get_derivative3_value(double phi) const;
    double get_derivative4_value(double phi) const;
    std::set<double> get_minimum_argument() const;
    std::set<double> get_minimum_argument_numerical() const;
    utility::Result<std::set<double>, NoKnownAnalyticalSolutionError> get_minimum_argument_analytical() const;

   private:
    static double calculate_prefactor(double cos_coef, double sin_coef, double sqcos_coef, double free_coef);
    AcosPlusBsinPlusCsqcosPlusZ(double cos_coef, double sin_coef, double sqcos_coef, double free_coef);
    const double _prefactor;
    const double _cos_coef;
    const double _sin_coef;
    const double _sqcos_coef;
    const double _free_coef;
    double get_value_without_prefactor(double phi) const;
    double get_derivative_value_without_prefactor(double phi) const;
    double get_derivative2_value_without_prefactor(double phi) const;
    double get_derivative3_value_without_prefactor(double phi) const;
    double get_derivative4_value_without_prefactor(double phi) const;
    bool is_degenerated_to_const_function() const;
    bool is_degenerated_to_sqcos_coef_equal_to_zero_case() const;
    bool is_degenerated_to_cos_coef_equal_to_zero_case() const;
    bool is_degenerated_to_sin_coef_equal_to_zero_case() const;
    std::set<double> get_minimum_argument_analytical_when_sqcos_coef_is_zero() const;
    std::set<double> get_minimum_argument_analytical_when_sqcos_coef_is_not_zero_and_cos_coef_is_zero() const;
    std::set<double> get_minimum_argument_analytical_when_sqcos_coef_is_not_zero_and_sin_coef_is_zero() const;
};

}  // end of namespace monostar_hamiltonians
