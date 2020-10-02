#include<model_monostar/hamiltonian_fo_params_helpers.hpp>

#include<model_monostar/simple_numerical_function_analyzer.hpp>

#include<utility/almost_equal.hpp>

#include<cmath>
#include<cassert>

#include<iostream> //TODO remove (debug)

// #######################################################################
// ## NoKnownAnalicycalSolutionError                                    ##
// #######################################################################

NoKnownAnalicycalSolutionError::NoKnownAnalicycalSolutionError() :
    std::domain_error("The params are out of the scope of analitical minimum finder.") {
}

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

// #######################################################################

BsinPlusCsqcosPlusZ BsinPlusCsqcosPlusZ::Builder::build() const {
    return BsinPlusCsqcosPlusZ(_sin_coef, _sqcos_coef, _free_coef);
}

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
        return get_minimum_argument_analitycal_when_sqcos_coef_is_zero();
    } else {
        return get_minimum_argument_analitycal_when_sqcos_coef_is_not_zero();
    }
}

std::set<double> BsinPlusCsqcosPlusZ::get_minimum_argument_analitycal_when_sqcos_coef_is_zero() const {
    assert(utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    if (_sin_coef < 0) {
        return {+M_PI/2};
    } else if (_sin_coef > 0) {
        return {-M_PI/2};
    } else {
        return {};
    }
}

std::set<double> BsinPlusCsqcosPlusZ::get_minimum_argument_analitycal_when_sqcos_coef_is_not_zero() const {
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
            std::cout << "HOOK AD" << std::endl; //TODO remove
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

// #######################################################################
// ## AcosPlusBsinPlusCsqcosPlusZ                                       ##
// #######################################################################

double AcosPlusBsinPlusCsqcosPlusZ::calculate_prefactor(double cos_coef, double sin_coef, double sqcos_coef, double free_coef) {
    double prefactor = std::sqrt(
                cos_coef * cos_coef +
                sin_coef * sin_coef +
                sqcos_coef * sqcos_coef +
                free_coef * free_coef);
    return (prefactor > 100 * std::numeric_limits<double>::epsilon() ?
                prefactor :
                1.0);
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

AcosPlusBsinPlusCsqcosPlusZ::AcosPlusBsinPlusCsqcosPlusZ(double cos_coef, double sin_coef, double sqcos_coef, double free_coef) :
    _prefactor(calculate_prefactor(cos_coef, sin_coef, sqcos_coef, free_coef)),
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
    return + _cos_coef * std::cos(phi)
            + _sin_coef * std::sin(phi)
            + _sqcos_coef * std::cos(phi) * std::cos(phi)
            + _free_coef;
}

double AcosPlusBsinPlusCsqcosPlusZ::get_derivative_value_without_prefactor(double phi) const {
    return - _cos_coef * std::sin(phi)
            + _sin_coef * std::cos(phi)
            - 2 * _sqcos_coef * std::cos(phi) * std::sin(phi);
}

double AcosPlusBsinPlusCsqcosPlusZ::get_value(double phi) const {
    return _prefactor * get_value_without_prefactor(phi);
}

double AcosPlusBsinPlusCsqcosPlusZ::get_derivative_value(double phi) const {
    return _prefactor * get_derivative_value_without_prefactor(phi);
}

std::set<double> AcosPlusBsinPlusCsqcosPlusZ::get_minimum_argument() const {
    if (const auto& analicical_solution = get_minimum_argument_analitycal()) {
        return analicical_solution.unwrap();
    }
    const auto& numerical_solution = get_minimum_argument_numerical();
    return numerical_solution;
}

std::set<double> AcosPlusBsinPlusCsqcosPlusZ::get_minimum_argument_analitycal_when_sqcos_coef_is_zero() const {
    assert(utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    return AcosPlusBsinPlusZ::Builder()
            .set_cos_coef(_cos_coef)
            .set_sin_coef(_sin_coef)
            .set_free_coef(_free_coef)
            .build()
            .get_minimum_argument();
}

std::set<double> AcosPlusBsinPlusCsqcosPlusZ::get_minimum_argument_analitycal_when_sqcos_coef_is_not_zero_and_cos_coef_is_zero() const {
    assert(!utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    assert(utility::almost_equal(_cos_coef, 0.0, 1.0));
    return BsinPlusCsqcosPlusZ::Builder()
            .set_sin_coef(_sin_coef)
            .set_sqcos_coef(_sqcos_coef)
            .set_free_coef(_free_coef)
            .build()
            .get_minimum_argument();
}

std::set<double> AcosPlusBsinPlusCsqcosPlusZ::get_minimum_argument_analitycal_when_sqcos_coef_is_not_zero_and_sin_coef_is_zero() const {
    assert(!utility::almost_equal(_sqcos_coef, 0.0, 1.0));
    assert(utility::almost_equal(_sin_coef, 0.0, 1.0));
    assert(false);
    //TODO implement
    return {0.0/0.0};
}

std::set<double> AcosPlusBsinPlusCsqcosPlusZ::get_minimum_argument_numerical() const {
    const std::function<double(double)> fn = [this](double phi){return get_value_without_prefactor(phi);};
    const std::function<double(double)> fn_prim = [this](double phi){return get_derivative_value_without_prefactor(phi);};
    return find_all_global_minima_periodic_2_pi(fn, fn_prim);
}

utility::Result<std::set<double>, NoKnownAnalicycalSolutionError> AcosPlusBsinPlusCsqcosPlusZ::get_minimum_argument_analitycal() const {
    using ResultT = utility::Result<std::set<double>, NoKnownAnalicycalSolutionError>;
    if (is_degenerated_to_const_function()) {
        return ResultT::Ok({});
    } else if (is_degenerated_to_sqcos_coef_equal_to_zero_case()) {
        return ResultT::Ok(get_minimum_argument_analitycal_when_sqcos_coef_is_zero());
    } else if (is_degenerated_to_cos_coef_equal_to_zero_case()) {
        return ResultT::Ok(get_minimum_argument_analitycal_when_sqcos_coef_is_not_zero_and_cos_coef_is_zero());
    } else if (is_degenerated_to_sin_coef_equal_to_zero_case()) {
        return ResultT::Ok(get_minimum_argument_analitycal_when_sqcos_coef_is_not_zero_and_sin_coef_is_zero());
    } else {
        return ResultT::Err(NoKnownAnalicycalSolutionError());
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
