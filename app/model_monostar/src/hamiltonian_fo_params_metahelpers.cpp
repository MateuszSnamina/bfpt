#include<model_monostar/hamiltonian_fo_params_metahelpers.hpp>

#include<model_monostar/simple_numerical_function_analyzer.hpp>

#include<utility/almost_equal.hpp>

#include<cmath>
#include<cassert>

#include<iostream> //TODO remove (debug)

// #######################################################################
// ## AcosPlusBsinPlusZ                                                 ##
// #######################################################################

AcosPlusBsinPlusZ::Builder AcosPlusBsinPlusZ::Builder::set_cos_coef(double cos_coef) {
    _cos_coef = cos_coef;
    return *this;
}

AcosPlusBsinPlusZ::Builder AcosPlusBsinPlusZ::Builder::set_sin_coef(double sin_coef) {
    _sin_coef = sin_coef;
    return *this;
}

AcosPlusBsinPlusZ::Builder AcosPlusBsinPlusZ::Builder::set_free_coef(double free_coef) {
    _free_coef = free_coef;
    return *this;
}

AcosPlusBsinPlusZ AcosPlusBsinPlusZ::Builder::build() const {
    return AcosPlusBsinPlusZ(_cos_coef, _sin_coef, _free_coef);
}

// #######################################################################

AcosPlusBsinPlusZ::AcosPlusBsinPlusZ(double cos_coef, double sin_coef, double free_coef) :
    _cos_coef(cos_coef),
    _sin_coef(sin_coef),
    _free_coef(free_coef) {
}

double AcosPlusBsinPlusZ::get_cos_coef() const {
    return _cos_coef;
}

double AcosPlusBsinPlusZ::get_sin_coef() const {
    return _sin_coef;
}

double AcosPlusBsinPlusZ::get_free_coef() const {
    return _free_coef;
}

double AcosPlusBsinPlusZ::get_value(double phi) const {
    return + _cos_coef * std::cos(phi)
            + _sin_coef * std::sin(phi)
            + _free_coef;
}

std::set<double> AcosPlusBsinPlusZ::get_minimum_argument() const {
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

bool AcosPlusBsinPlusZ::is_degenerated_to_const_function() const {
    const double h = std::hypot(_cos_coef, _sin_coef);
    return h < 100 * std::numeric_limits<double>::epsilon();
}

// #######################################################################
// ## BsinPlusCsqcosPlusZ                                               ##
// #######################################################################

BsinPlusCsqcosPlusZ::Builder BsinPlusCsqcosPlusZ::Builder::set_sin_coef(double sin_coef) {
    _sin_coef = sin_coef;
    return *this;
}

BsinPlusCsqcosPlusZ::Builder BsinPlusCsqcosPlusZ::Builder::set_sqcos_coef(double sqcos_coef) {
    _sqcos_coef = sqcos_coef;
    return *this;
}

BsinPlusCsqcosPlusZ::Builder BsinPlusCsqcosPlusZ::Builder::set_free_coef(double free_coef) {
    _free_coef = free_coef;
    return *this;
}

BsinPlusCsqcosPlusZ BsinPlusCsqcosPlusZ::Builder::build() const {
    return BsinPlusCsqcosPlusZ(_sin_coef, _sqcos_coef, _free_coef);
}

// #######################################################################

BsinPlusCsqcosPlusZ::BsinPlusCsqcosPlusZ(double sin_coef, double sqcos_coef, double free_coef) :
    _sin_coef(sin_coef),
    _sqcos_coef(sqcos_coef),
    _free_coef(free_coef) {
}


double BsinPlusCsqcosPlusZ::get_sin_coef() const {
    return _sin_coef;
}

double BsinPlusCsqcosPlusZ::get_sqcos_coef() const {
    return _sqcos_coef;
}

double BsinPlusCsqcosPlusZ::get_free_coef() const {
    return _free_coef;
}

double BsinPlusCsqcosPlusZ::get_value(double phi) const {
    return + _sin_coef * std::sin(phi)
            + _sqcos_coef * std::cos(phi) * std::cos(phi)
            + _free_coef;
}

std::set<double> BsinPlusCsqcosPlusZ::get_minimum_argument() const {
    if (is_degenerated_to_const_function()) {
        return {};
    } else if (is_degenerated_to_sqcos_coef_equal_to_zero_case()) {
        return get_minimum_argument_analytical_when_sqcos_coef_is_zero();
    } else {
        return get_minimum_argument_analytical_when_sqcos_coef_is_not_zero();
    }
}

std::set<double> BsinPlusCsqcosPlusZ::get_minimum_argument_analytical_when_sqcos_coef_is_zero() const {
    assert(utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    if (_sin_coef < 0) {
        return {+M_PI/2};
    } else if (_sin_coef > 0) {
        return {-M_PI/2};
    } else {
        return {};
    }
}

std::set<double> BsinPlusCsqcosPlusZ::get_minimum_argument_analytical_when_sqcos_coef_is_not_zero() const {
    assert(!utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    if (_sin_coef <= -std::abs(2 * _sqcos_coef)) {
        return {+M_PI/2};
    } else if (_sin_coef >= +std::abs(2 * _sqcos_coef)) {
        return {-M_PI/2};
    } else {
        if (_sqcos_coef > 0) {
            assert(std::abs(_sin_coef) < std::abs(2 * _sqcos_coef));
            const double theta_1 = +M_PI/2;
            const double theta_2 = -M_PI/2;
            [[maybe_unused]] const double f_theta_1 = get_value(theta_1);
            [[maybe_unused]] const double f_theta_2 = get_value(theta_2);
            assert(utility::almost_equal(f_theta_1, +_sin_coef + _free_coef, 1.0));
            assert(utility::almost_equal(f_theta_2, -_sin_coef + _free_coef, 1.0));
            if (_sin_coef < 0) {
                assert(f_theta_1 < f_theta_2 || utility::almost_equal(f_theta_1, f_theta_2, 1.0));
                return {theta_1};
            } else if (_sin_coef > 0) {
                assert(f_theta_2 < f_theta_1 || utility::almost_equal(f_theta_1, f_theta_2, 1.0));
                return {theta_2};
            } else {
                assert(_sin_coef == 0);
                assert(utility::almost_equal(f_theta_1, f_theta_2, 1.0));
                return {theta_1, theta_2};
            }
        } else {
            assert(std::abs(_sin_coef) < std::abs(2 * _sqcos_coef));
            assert(_sqcos_coef < 0);
            const double B_prim = _sin_coef/ (2 * _sqcos_coef);
            const double theta_3 = std::asin(B_prim);
            const double theta_4_pre = M_PI - std::asin(B_prim);
            const double theta_4 = (theta_4_pre <= M_PI ? theta_4_pre : theta_4_pre - 2 * M_PI);
            [[maybe_unused]] const double f_theta_3 = get_value(theta_3);
            [[maybe_unused]] const double f_theta_4 = get_value(theta_4);
            assert(utility::almost_equal(f_theta_3, f_theta_4, 1.0));
            return {theta_3, theta_4};
        }
    }
}

bool BsinPlusCsqcosPlusZ::is_degenerated_to_const_function() const {
    const double h = std::sqrt(_sqcos_coef * _sqcos_coef + _sin_coef * _sin_coef);
    return h < 100 * std::numeric_limits<double>::epsilon();
}

bool BsinPlusCsqcosPlusZ::is_degenerated_to_sqcos_coef_equal_to_zero_case() const {
    const double N = std::abs(_sin_coef);
    return std::abs(_sqcos_coef) < 100 * N * std::numeric_limits<double>::epsilon();
}

//// #######################################################################
//// ## AcosPlusCsqcosPlusZ                                               ##
//// #######################################################################

AcosPlusCsqcosPlusZ::Builder AcosPlusCsqcosPlusZ::Builder::set_cos_coef(double cos_coef) {
    _cos_coef = cos_coef;
    return *this;
}

AcosPlusCsqcosPlusZ::Builder AcosPlusCsqcosPlusZ::Builder::set_sqcos_coef(double sqcos_coef) {
    _sqcos_coef = sqcos_coef;
    return *this;
}

AcosPlusCsqcosPlusZ::Builder AcosPlusCsqcosPlusZ::Builder::set_free_coef(double free_coef) {
    _free_coef = free_coef;
    return *this;
}

AcosPlusCsqcosPlusZ AcosPlusCsqcosPlusZ::Builder::build() const {
    return AcosPlusCsqcosPlusZ(_cos_coef, _sqcos_coef, _free_coef);
}

// #######################################################################

AcosPlusCsqcosPlusZ::AcosPlusCsqcosPlusZ(double cos_coef, double sqcos_coef, double free_coef) :
    _cos_coef(cos_coef),
    _sqcos_coef(sqcos_coef),
    _free_coef(free_coef) {
}


double AcosPlusCsqcosPlusZ::get_cos_coef() const {
    return _cos_coef;
}

double AcosPlusCsqcosPlusZ::get_sqcos_coef() const {
    return _sqcos_coef;
}

double AcosPlusCsqcosPlusZ::get_free_coef() const {
    return _free_coef;
}

double AcosPlusCsqcosPlusZ::get_value(double phi) const {
    return + _cos_coef * std::cos(phi)
            + _sqcos_coef * std::cos(phi) * std::cos(phi)
            + _free_coef;
}

std::set<double> AcosPlusCsqcosPlusZ::get_minimum_argument() const {
    if (is_degenerated_to_const_function()) {
        return {};
    } else if (is_degenerated_to_sqcos_coef_equal_to_zero_case()) {
        return get_minimum_argument_analytical_when_sqcos_coef_is_zero();
    } else {
        return get_minimum_argument_analytical_when_sqcos_coef_is_not_zero();
    }
}

std::set<double> AcosPlusCsqcosPlusZ::get_minimum_argument_analytical_when_sqcos_coef_is_zero() const {
    assert(utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    if (_cos_coef < 0) {
        return {+0};
    } else if (_cos_coef > 0) {
        return {+M_PI};
    } else {
        return {};
    }
}

std::set<double> AcosPlusCsqcosPlusZ::get_minimum_argument_analytical_when_sqcos_coef_is_not_zero() const {
    assert(!utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    if (_cos_coef <= -std::abs(2 * _sqcos_coef)) {
        return {0};
    } else if (_cos_coef >= +std::abs(2 * _sqcos_coef)) {
        return {M_PI};
    } else {
        if (_sqcos_coef < 0) {
            assert(std::abs(_cos_coef) < std::abs(2 * _sqcos_coef));
            const double theta_1 = 0;
            const double theta_2 = M_PI;
            [[maybe_unused]] const double f_theta_1 = get_value(theta_1);
            [[maybe_unused]] const double f_theta_2 = get_value(theta_2);
            assert(utility::almost_equal(f_theta_1, +_cos_coef + _sqcos_coef + _free_coef, 1.0));
            assert(utility::almost_equal(f_theta_2, -_cos_coef + _sqcos_coef + _free_coef, 1.0));
            if (_cos_coef < 0) {
                assert(f_theta_1 < f_theta_2 || utility::almost_equal(f_theta_1, f_theta_2, 1.0));
                return {theta_1};
            } else if (_cos_coef > 0) {
                assert(f_theta_2 < f_theta_1 || utility::almost_equal(f_theta_1, f_theta_2, 1.0));
                return {theta_2};
            } else {
                assert(_cos_coef == 0);
                assert(utility::almost_equal(f_theta_1, f_theta_2, 1.0));
                return {theta_1, theta_2};
            }
        } else {
            assert(std::abs(_cos_coef) < std::abs(2 * _sqcos_coef));
            assert(_sqcos_coef > 0);
            const double A_prim = -_cos_coef/ (2 * _sqcos_coef);
            const double theta_3 = + std::acos(A_prim);
            const double theta_4 = - std::acos(A_prim);
            [[maybe_unused]] const double f_theta_3 = get_value(theta_3);
            [[maybe_unused]] const double f_theta_4 = get_value(theta_4);
            assert(utility::almost_equal(f_theta_3, f_theta_4, 1.0));
            return {theta_3, theta_4};
        }
    }
}

bool AcosPlusCsqcosPlusZ::is_degenerated_to_const_function() const {
    const double h = std::sqrt(_sqcos_coef * _sqcos_coef + _cos_coef * _cos_coef);
    return h < 100 * std::numeric_limits<double>::epsilon();
}

bool AcosPlusCsqcosPlusZ::is_degenerated_to_sqcos_coef_equal_to_zero_case() const {
    const double N = std::abs(_cos_coef);
    return std::abs(_sqcos_coef) < 100 * N * std::numeric_limits<double>::epsilon();
}
