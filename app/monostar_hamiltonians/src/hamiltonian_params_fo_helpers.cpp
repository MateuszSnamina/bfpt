#include <monostar_hamiltonians/hamiltonian_params_fo_helpers.hpp>

#include <monostar_hamiltonians/hamiltonian_params_fo_metahelpers.hpp>
#include <monostar_hamiltonians/simple_numerical_function_analyzer.hpp>

#include <utility/almost_equal.hpp>

#include <cmath>
#include <cassert>

// #######################################################################
// ## NoKnownAnalyticalSolutionError                                    ##
// #######################################################################

namespace monostar_hamiltonians {

NoKnownAnalyticalSolutionError::NoKnownAnalyticalSolutionError() : std::domain_error("The params are out of the scope of analytical minimum finder.") {
}

}  // end of namespace monostar_hamiltonians

// #######################################################################
// ## AcosPlusBsinPlusCsqcosPlusZ                                       ##
// #######################################################################

namespace monostar_hamiltonians {

double AcosPlusBsinPlusCsqcosPlusZ::calculate_prefactor(double cos_coef, double sin_coef, double sqcos_coef, double free_coef) {
    double prefactor = std::sqrt(
        cos_coef * cos_coef +
        sin_coef * sin_coef +
        sqcos_coef * sqcos_coef +
        free_coef * free_coef);
    return (prefactor > 100 * std::numeric_limits<double>::epsilon() ? prefactor : 1.0);
}

AcosPlusBsinPlusCsqcosPlusZ::Builder AcosPlusBsinPlusCsqcosPlusZ::Builder::set_cos_coef(double cos_coef) {
    _cos_coef = cos_coef;
    return *this;
}

AcosPlusBsinPlusCsqcosPlusZ::Builder AcosPlusBsinPlusCsqcosPlusZ::Builder::set_sin_coef(double sin_coef) {
    _sin_coef = sin_coef;
    return *this;
}

AcosPlusBsinPlusCsqcosPlusZ::Builder AcosPlusBsinPlusCsqcosPlusZ::Builder::set_sqcos_coef(double sqcos_coef) {
    _sqcos_coef = sqcos_coef;
    return *this;
}

AcosPlusBsinPlusCsqcosPlusZ::Builder AcosPlusBsinPlusCsqcosPlusZ::Builder::set_free_coef(double free_coef) {
    _free_coef = free_coef;
    return *this;
}

// #######################################################################

AcosPlusBsinPlusCsqcosPlusZ AcosPlusBsinPlusCsqcosPlusZ::Builder::build() const {
    return AcosPlusBsinPlusCsqcosPlusZ(_cos_coef, _sin_coef, _sqcos_coef, _free_coef);
}

AcosPlusBsinPlusCsqcosPlusZ::AcosPlusBsinPlusCsqcosPlusZ(double cos_coef, double sin_coef, double sqcos_coef, double free_coef)
    : _prefactor(calculate_prefactor(cos_coef, sin_coef, sqcos_coef, free_coef)),
      _cos_coef(cos_coef / _prefactor),
      _sin_coef(sin_coef / _prefactor),
      _sqcos_coef(sqcos_coef / _prefactor),
      _free_coef(free_coef / _prefactor) {
}

double AcosPlusBsinPlusCsqcosPlusZ::get_cos_coef() const {
    return _prefactor * _cos_coef;
}

double AcosPlusBsinPlusCsqcosPlusZ::get_sin_coef() const {
    return _prefactor * _sin_coef;
}

double AcosPlusBsinPlusCsqcosPlusZ::get_sqcos_coef() const {
    return _prefactor * _sqcos_coef;
}

double AcosPlusBsinPlusCsqcosPlusZ::get_free_coef() const {
    return _prefactor * _free_coef;
}

double AcosPlusBsinPlusCsqcosPlusZ::get_value_without_prefactor(double phi) const {
    return +_cos_coef * std::cos(phi) + _sin_coef * std::sin(phi) + _sqcos_coef * std::cos(phi) * std::cos(phi) + _free_coef;
}

double AcosPlusBsinPlusCsqcosPlusZ::get_derivative_value_without_prefactor(double phi) const {
    return -_cos_coef * std::sin(phi) + _sin_coef * std::cos(phi) - 2 * _sqcos_coef * std::cos(phi) * std::sin(phi);
}

double AcosPlusBsinPlusCsqcosPlusZ::get_derivative2_value_without_prefactor(double phi) const {
    return -_cos_coef * std::cos(phi) - _sin_coef * std::sin(phi) - 4 * _sqcos_coef * std::cos(phi) * std::cos(phi) + 2 * _sqcos_coef;
}

double AcosPlusBsinPlusCsqcosPlusZ::get_derivative3_value_without_prefactor(double phi) const {
    return +_cos_coef * std::sin(phi) - _sin_coef * std::cos(phi) + 8 * _sqcos_coef * std::cos(phi) * std::sin(phi);
}

double AcosPlusBsinPlusCsqcosPlusZ::get_derivative4_value_without_prefactor(double phi) const {
    return +_cos_coef * std::cos(phi) + _sin_coef * std::sin(phi) + 16 * _sqcos_coef * std::cos(phi) * std::cos(phi) - 8 * _sqcos_coef;
}

double AcosPlusBsinPlusCsqcosPlusZ::get_value(double phi) const {
    return _prefactor * get_value_without_prefactor(phi);
}

double AcosPlusBsinPlusCsqcosPlusZ::get_derivative_value(double phi) const {
    return _prefactor * get_derivative_value_without_prefactor(phi);
}

double AcosPlusBsinPlusCsqcosPlusZ::get_derivative2_value(double phi) const {
    return _prefactor * get_derivative2_value_without_prefactor(phi);
}

double AcosPlusBsinPlusCsqcosPlusZ::get_derivative3_value(double phi) const {
    return _prefactor * get_derivative3_value_without_prefactor(phi);
}

double AcosPlusBsinPlusCsqcosPlusZ::get_derivative4_value(double phi) const {
    return _prefactor * get_derivative4_value_without_prefactor(phi);
}

std::set<double> AcosPlusBsinPlusCsqcosPlusZ::get_minimum_argument() const {
    if (const auto& analicical_solution = get_minimum_argument_analytical()) {
        return analicical_solution.unwrap();
    }
    const auto& numerical_solution = get_minimum_argument_numerical();
    return numerical_solution;
}

std::set<double> AcosPlusBsinPlusCsqcosPlusZ::get_minimum_argument_analytical_when_sqcos_coef_is_zero() const {
    assert(utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    return AcosPlusBsinPlusZ::Builder()
        .set_cos_coef(_cos_coef)
        .set_sin_coef(_sin_coef)
        .set_free_coef(_free_coef)
        .build()
        .get_minimum_argument();
}

std::set<double> AcosPlusBsinPlusCsqcosPlusZ::get_minimum_argument_analytical_when_sqcos_coef_is_not_zero_and_cos_coef_is_zero() const {
    assert(!utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    assert(utility::almost_equal(_cos_coef, 0.0, 1.0));
    return BsinPlusCsqcosPlusZ::Builder()
        .set_sin_coef(_sin_coef)
        .set_sqcos_coef(_sqcos_coef)
        .set_free_coef(_free_coef)
        .build()
        .get_minimum_argument();
}

std::set<double> AcosPlusBsinPlusCsqcosPlusZ::get_minimum_argument_analytical_when_sqcos_coef_is_not_zero_and_sin_coef_is_zero() const {
    assert(!utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    assert(utility::almost_equal(_sin_coef, 0.0, 1.0));
    return AcosPlusCsqcosPlusZ::Builder()
        .set_cos_coef(_cos_coef)
        .set_sqcos_coef(_sqcos_coef)
        .set_free_coef(_free_coef)
        .build()
        .get_minimum_argument();
}

std::set<double> AcosPlusBsinPlusCsqcosPlusZ::get_minimum_argument_numerical() const {
    const std::function<double(double)> fn = [this](double phi) { return get_value_without_prefactor(phi); };
    const std::function<double(double)> fn_prim = [this](double phi) { return get_derivative_value_without_prefactor(phi); };
    return function_analyzer::find_all_global_minima_periodic_2_pi(fn, fn_prim);
}

utility::Result<std::set<double>, NoKnownAnalyticalSolutionError> AcosPlusBsinPlusCsqcosPlusZ::get_minimum_argument_analytical() const {
    using ResultT = utility::Result<std::set<double>, NoKnownAnalyticalSolutionError>;
    if (is_degenerated_to_const_function()) {
        return ResultT::Ok({});
    } else if (is_degenerated_to_sqcos_coef_equal_to_zero_case()) {
        return ResultT::Ok(get_minimum_argument_analytical_when_sqcos_coef_is_zero());
    } else if (is_degenerated_to_cos_coef_equal_to_zero_case()) {
        return ResultT::Ok(get_minimum_argument_analytical_when_sqcos_coef_is_not_zero_and_cos_coef_is_zero());
    } else if (is_degenerated_to_sin_coef_equal_to_zero_case()) {
        return ResultT::Ok(get_minimum_argument_analytical_when_sqcos_coef_is_not_zero_and_sin_coef_is_zero());
    } else {
        return ResultT::Err(NoKnownAnalyticalSolutionError());
    }
}

bool AcosPlusBsinPlusCsqcosPlusZ::is_degenerated_to_const_function() const {
    const double h = std::sqrt(_cos_coef * _cos_coef + _sqcos_coef * _sqcos_coef + _sin_coef * _sin_coef);
    return h < 100 * std::numeric_limits<double>::epsilon();
}

bool AcosPlusBsinPlusCsqcosPlusZ::is_degenerated_to_cos_coef_equal_to_zero_case() const {
    const double N = std::hypot(_sqcos_coef, _sin_coef);
    return std::abs(_cos_coef) < 100 * N * std::numeric_limits<double>::epsilon();
}

bool AcosPlusBsinPlusCsqcosPlusZ::is_degenerated_to_sin_coef_equal_to_zero_case() const {
    const double N = std::hypot(_cos_coef, _sqcos_coef);
    return std::abs(_sin_coef) < 100 * N * std::numeric_limits<double>::epsilon();
}

bool AcosPlusBsinPlusCsqcosPlusZ::is_degenerated_to_sqcos_coef_equal_to_zero_case() const {
    const double N = std::hypot(_cos_coef, _sin_coef);
    return std::abs(_sqcos_coef) < 100 * N * std::numeric_limits<double>::epsilon();
}

// #######################################################################

AcosPlusBsinPlusCsqcosPlusZ operator+(const AcosPlusBsinPlusCsqcosPlusZ& fun1, const AcosPlusBsinPlusCsqcosPlusZ& fun2) {
    return AcosPlusBsinPlusCsqcosPlusZ::Builder()
        .set_cos_coef(fun1.get_cos_coef() + fun2.get_cos_coef())
        .set_sin_coef(fun1.get_sin_coef() + fun2.get_sin_coef())
        .set_sqcos_coef(fun1.get_sqcos_coef() + fun2.get_sqcos_coef())
        .set_free_coef(fun1.get_free_coef() + fun2.get_free_coef())
        .build();
}

AcosPlusBsinPlusCsqcosPlusZ operator*(double factor, const AcosPlusBsinPlusCsqcosPlusZ& fun) {
    return AcosPlusBsinPlusCsqcosPlusZ::Builder()
        .set_cos_coef(factor * fun.get_cos_coef())
        .set_sin_coef(factor * fun.get_sin_coef())
        .set_sqcos_coef(factor * fun.get_sqcos_coef())
        .set_free_coef(factor * fun.get_free_coef())
        .build();
}

}  // end of namespace monostar_hamiltonians
