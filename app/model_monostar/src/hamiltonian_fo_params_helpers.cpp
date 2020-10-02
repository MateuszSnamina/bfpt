#include<model_monostar/hamiltonian_fo_params_helpers.hpp>

#include<model_monostar/simple_numerical_function_analyzer.hpp>

#include<utility/almost_equal.hpp>

#include<cmath>
#include<cassert>

#include<iostream> //TODO remove (debug)

// #######################################################################
// ## AcosPlucBsinPlusZ                                                 ##
// #######################################################################

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

// #######################################################################

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

std::set<double> AcosPlucBsinPlusZ::get_minimum_argument() const {
    const double h = std::hypot(_cos_coef, _sin_coef);
    if (is_degenerated_to_const_function()) {
        return {};
    }
    const double x = - _cos_coef / h;
    const double y = - _sin_coef / h;
    const double phi_0 = std::atan2(y, x);
    if (std::abs(phi_0 + M_PI) < std::numeric_limits<double>::epsilon()) {
        return {M_PI};
    } else {
        return {phi_0};
    }
}

bool AcosPlucBsinPlusZ::is_degenerated_to_const_function() const {
    const double h = std::hypot(_cos_coef, _sin_coef);
    return h < 100 * std::numeric_limits<double>::epsilon();
}

// #######################################################################
// ## AcosPlucBsinPlusCsqcosPlusZ                                       ##
// #######################################################################

double AcosPlucBsinPlusCsqcosPlusZ::calculate_prefactor(double cos_coef, double sin_coef, double sqcos_coef, double free_coef) {
    double prefactor = std::sqrt(
                cos_coef * cos_coef +
                sin_coef * sin_coef +
                sqcos_coef * sqcos_coef +
                free_coef * free_coef);
    return (prefactor > 100 * std::numeric_limits<double>::epsilon() ?
                prefactor :
                1.0);
}

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

// #######################################################################

AcosPlucBsinPlusCsqcosPlusZ AcosPlucBsinPlusCsqcosPlusZ::Builder::build() const {
    return AcosPlucBsinPlusCsqcosPlusZ(_cos_coef, _sin_coef, _sqcos_coef, _free_coef);
}

AcosPlucBsinPlusCsqcosPlusZ::AcosPlucBsinPlusCsqcosPlusZ(double cos_coef, double sin_coef, double sqcos_coef, double free_coef) :
    _prefactor(calculate_prefactor(cos_coef, sin_coef, sqcos_coef, free_coef)),
    _cos_coef(cos_coef / _prefactor),
    _sin_coef(sin_coef / _prefactor),
    _sqcos_coef(sqcos_coef / _prefactor),
    _free_coef(free_coef / _prefactor) {
}

double AcosPlucBsinPlusCsqcosPlusZ::get_cos_coef() const {
    return _prefactor * _cos_coef;
}

double AcosPlucBsinPlusCsqcosPlusZ::get_sin_coef() const {
    return _prefactor * _sin_coef;
}

double AcosPlucBsinPlusCsqcosPlusZ::get_sqcos_coef() const {
    return _prefactor * _sqcos_coef;
}

double AcosPlucBsinPlusCsqcosPlusZ::get_free_coef() const {
    return _prefactor * _free_coef;
}


double AcosPlucBsinPlusCsqcosPlusZ::get_value_without_prefactor(double phi) const {
    return + _cos_coef * std::cos(phi)
            + _sin_coef * std::sin(phi)
            + _sqcos_coef * std::cos(phi) * std::cos(phi)
            + _free_coef;
}

double AcosPlucBsinPlusCsqcosPlusZ::get_derivative_value_without_prefactor(double phi) const {
    return - _cos_coef * std::sin(phi)
            + _sin_coef * std::cos(phi)
            - 2 * _sqcos_coef * std::cos(phi) * std::sin(phi);
}

double AcosPlucBsinPlusCsqcosPlusZ::get_value(double phi) const {
    return _prefactor * get_value_without_prefactor(phi);
}

double AcosPlucBsinPlusCsqcosPlusZ::get_derivative_value(double phi) const {
    return _prefactor * get_derivative_value_without_prefactor(phi);
}

std::set<double> AcosPlucBsinPlusCsqcosPlusZ::get_minimum_argument() const {
    { //TODO: remove
        std::cout << std::endl;
        std::cout << "[AcosPlucBsinPlusCsqcosPlusZ::get_minimum_argument] [Debug]: get_minimum_not_analitycal" << std::endl; //TODO: remove
        for  (const auto phi : get_minimum_not_analitycal()) { //TODO: remove
            std::cout << phi  << std::endl; //TODO: remove
        } //TODO: remove
        std::cout << "**************" << std::endl; //TODO: remove
    } //TODO: remove
    if (is_degenerated_to_const_function()) {
        return {};
    } else if (is_degenerated_to_sqcos_coef_equal_to_zero_case()) {
        { //TODO: remove
            std::cout << "[AcosPlucBsinPlusCsqcosPlusZ::get_minimum_argument] [Debug]: get_minimum_argument_when_cos_coef_is_zero" << std::endl; //TODO: remove
            for  (const auto phi : get_minimum_argument_when_sqcos_coef_is_zero()) { //TODO: remove
                std::cout << phi  << std::endl; //TODO: remove
            } //TODO: remove
            std::cout << "**************" << std::endl; //TODO: remove
        } //TODO: remove
        return get_minimum_argument_when_sqcos_coef_is_zero();
    } else if (is_degenerated_to_cos_coef_equal_to_zero_case()) {
        { //TODO: remove
            std::cout << "[AcosPlucBsinPlusCsqcosPlusZ::get_minimum_argument] [Debug]: get_minimum_argument_when_cos_coef_is_zero" << std::endl; //TODO: remove
            for  (const auto phi : get_minimum_argument_when_cos_coef_is_zero()) { //TODO: remove
                std::cout << phi  << std::endl; //TODO: remove
            } //TODO: remove
            std::cout << "**************" << std::endl; //TODO: remove
        } //TODO: remove
        return get_minimum_argument_when_cos_coef_is_zero();
    } else if (is_degenerated_to_sin_coef_equal_to_zero_case()) {
        return get_minimum_argument_when_sin_coef_is_zero();
    } else {
        return get_minimum_not_analitycal();
    }
}

std::set<double> AcosPlucBsinPlusCsqcosPlusZ::get_minimum_argument_when_cos_coef_is_zero() const {
    assert(!is_degenerated_to_const_function());
    assert(utility::almost_equal(_cos_coef, 0.0, 1.0));
    assert(!utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    if (_sin_coef <= -std::abs(2 * _sqcos_coef)) {
        return {+M_PI/2};
    } else if (_sin_coef >= +std::abs(2 * _sqcos_coef)) {
        return {-M_PI/2};
    } else {
        assert(std::abs(_sin_coef) < std::abs(2 * _sqcos_coef));
        if (_sqcos_coef > 0) {
            const double theta_1 = +M_PI/2;
            const double theta_2 = -M_PI/2;
            [[maybe_unused]] const double f_theta_1 = get_value_without_prefactor(theta_1);
            [[maybe_unused]] const double f_theta_2 = get_value_without_prefactor(theta_2);
            assert(utility::almost_equal(f_theta_1, +_sin_coef + _free_coef, 1.0));
            assert(utility::almost_equal(f_theta_2, -_sin_coef + _free_coef, 1.0));
            if (_sin_coef < 0) {
                assert(f_theta_1 < f_theta_2 || utility::almost_equal(f_theta_1, f_theta_2, 1.0));
                return {+theta_1};
            } else if (_sin_coef > 0) {
                assert(f_theta_2 < f_theta_1 || utility::almost_equal(f_theta_1, f_theta_2, 1.0));
                return {f_theta_2};
            } else {
                assert(_sin_coef == 0);
                assert(utility::almost_equal(f_theta_1, f_theta_2, 1.0));
                return {-M_PI/2, +M_PI/2};
            }
        } else {
            assert(std::abs(_sin_coef) < std::abs(2 * _sqcos_coef));
            assert(_sqcos_coef < 0);
            const double B_prim = _sin_coef/ (2 * _sqcos_coef);
            const double theta_3 = std::asin(B_prim);
            const double theta_4 = M_PI - std::asin(B_prim);
            [[maybe_unused]] const double f_theta_3 = get_value_without_prefactor(theta_3);
            [[maybe_unused]] const double f_theta_4 = get_value_without_prefactor(theta_4);
            assert(utility::almost_equal(f_theta_3, f_theta_4, 1.0));
            return {theta_3, theta_4};
        }
    }
}

std::set<double> AcosPlucBsinPlusCsqcosPlusZ::get_minimum_argument_when_sin_coef_is_zero() const {
    assert(!is_degenerated_to_const_function());
    assert(utility::almost_equal(_sin_coef, 0.0, 1.0));
    assert(!utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    assert(false);
    //TODO implement
    return {0.0/0.0};
}

std::set<double> AcosPlucBsinPlusCsqcosPlusZ::get_minimum_argument_when_sqcos_coef_is_zero() const {
    return AcosPlucBsinPlusZ::Builder()
            .set_cos_coef(_cos_coef)
            .set_sin_coef(_sin_coef)
            .set_free_coef(_free_coef)
            .build()
            .get_minimum_argument();
}

std::set<double> AcosPlucBsinPlusCsqcosPlusZ::get_minimum_not_analitycal() const {
    const std::function<double(double)> fn = [this](double phi){return get_value_without_prefactor(phi);};
    const std::function<double(double)> fn_prim = [this](double phi){return get_derivative_value_without_prefactor(phi);};
    return find_all_global_minima_periodic_2_pi(fn, fn_prim);
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
