#include<monostar_hamiltonians/hamiltonian_params_fo_site_matrices.hpp>

#include<cmath>
#include<cassert>


// #######################################################################
// ## OneSiteOrbitalMatrices                                            ##
// #######################################################################

namespace monostar_hamiltonians {

arma::cx_mat22 OneSiteOrbitalMatrices::get_P_z_in_zx_basis() {
    const arma::mat22 P_z_in_zx_basis_re = {
        {1.0, 0.0},
        {0.0, 0.0}};
    arma::cx_mat P_z_in_zx_basis{P_z_in_zx_basis_re, arma::zeros(2, 2)};
    return P_z_in_zx_basis;
}


arma::cx_mat22 OneSiteOrbitalMatrices::get_P_x_in_zx_basis() {
    const arma::mat22 P_x_in_zx_basis_re = {
        {0.0, 0.0},
        {0.0, 1.0}};
    arma::cx_mat P_x_in_zx_basis{P_x_in_zx_basis_re, arma::zeros(2, 2)};
    return P_x_in_zx_basis;
}


arma::cx_mat22 OneSiteOrbitalMatrices::get_P_plus_in_zx_basis() {
    const arma::mat22 P_plus_in_zx_basis_re = {
        {+0.5, +0.5},
        {+0.5, +0.5}};
    arma::cx_mat P_plus_in_zx_basis{P_plus_in_zx_basis_re, arma::zeros(2, 2)};
    return P_plus_in_zx_basis;
}


arma::cx_mat22 OneSiteOrbitalMatrices::get_P_minus_in_zx_basis() {
    const arma::mat22 P_minus_in_zx_basis_re = {
        {+0.5, -0.5},
        {-0.5, +0.5}};
    arma::cx_mat P_minus_in_zx_basis{P_minus_in_zx_basis_re, arma::zeros(2, 2)};
    return P_minus_in_zx_basis;
}


arma::cx_mat22 OneSiteOrbitalMatrices::get_tau_z_in_zx_basis() {
    return OneSiteOrbitalMatrices::get_P_z_in_zx_basis() - OneSiteOrbitalMatrices::get_P_x_in_zx_basis();
}


arma::cx_mat22 OneSiteOrbitalMatrices::get_tau_x_in_zx_basis() {
    return OneSiteOrbitalMatrices::get_P_x_in_zx_basis() - OneSiteOrbitalMatrices::get_P_z_in_zx_basis();
}


arma::cx_mat22 OneSiteOrbitalMatrices::get_tau_plus_in_zx_basis() {
    return OneSiteOrbitalMatrices::get_P_plus_in_zx_basis() - OneSiteOrbitalMatrices::get_P_minus_in_zx_basis();
}


arma::cx_mat22 OneSiteOrbitalMatrices::get_tau_minus_in_zx_basis() {
    return OneSiteOrbitalMatrices::get_P_minus_in_zx_basis() - OneSiteOrbitalMatrices::get_P_plus_in_zx_basis();
}


/*
 * ╭         ╮     ╭         ╮              ╭                        ╮
 * │ |∥⟩ |⟂⟩ │  =  │ |x⟩ |z⟩ │ β, where β = │  +cos(θ/2)   -sin(θ/2) │
 * ╰         ╯     ╰         ╯              │  +sin(θ/2)   +cos(θ/2) │
 *                                          ╰                        ╯
 */
arma::cx_mat22 OneSiteOrbitalMatrices::get_beta_from_zx_to_ge(double orbital_theta) {
    const double s2 = std::sin(orbital_theta / 2.0), c2 = std::cos(orbital_theta / 2.0);
    const arma::mat22 beta_re = {
        {+c2, -s2},
        {+s2, +c2}
    };
    const arma::cx_mat22 beta{beta_re, arma::zeros(2, 2)};
    assert(arma::norm(beta.t() * beta - arma::eye(2, 2)) < 1e-4);
    assert(arma::norm(beta * beta.t() - arma::eye(2, 2)) < 1e-4);
    assert(std::abs(arma::det(beta) -1.0) < 1e-4);
    return beta;
}


/*
 *                 |∥⟩  |⟂⟩            |∥⟩      |⟂⟩
 *              ╭          ╮       ╭                 ╮
 * Pᶻ = ½ + ½ * │ +c1  -s1 │ |∥⟩ = │ +(c2)²   -s2*c2 │ |∥⟩
 *              │ -s1  -c1 │ |⟂⟩   │ -s2*c2   +(s2)² │ |⟂⟩
 *              ╰          ╯       ╰                 ╯
 */
arma::cx_mat22 OneSiteOrbitalMatrices::get_P_z_in_ge_basis(double orbital_theta) {
    const double s1 = std::sin(orbital_theta), c1 = std::cos(orbital_theta);
    const arma::mat22 P_z_in_ge_re_builder = {
        {+c1, -s1},
        {-s1, -c1}
    };
    const arma::mat22 P_z_in_ge_re = 0.5 * (arma::eye(2, 2) + P_z_in_ge_re_builder);
    const arma::cx_mat22 P_z_in_ge{P_z_in_ge_re, arma::zeros(2, 2)};
#ifndef NDEBUG
    const double s2 = std::sin(orbital_theta / 2.0), c2 = std::cos(orbital_theta / 2.0);
    const arma::mat22 P_z_in_ge_re_alternative{
        {+c2 * c2, -s2 * c2},
        {-s2 * c2, +s2 * s2}
    };
    assert(arma::norm(P_z_in_ge_re_alternative - P_z_in_ge_re) < 1e-4);
    const auto P_z_in_zx = OneSiteOrbitalMatrices::get_P_z_in_zx_basis();
    const auto beta = OneSiteOrbitalMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(P_z_in_ge - beta.t() * P_z_in_zx * beta) < 1e-4);
#endif
    return P_z_in_ge;
}


/*
 *                 |∥⟩  |⟂⟩            |∥⟩      |⟂⟩
 *              ╭          ╮       ╭                 ╮
 * Pˣ = ½ + ½ * │ -c1  +s1 │ |∥⟩ = │ +(s2)²   +s2*c2 │ |∥⟩
 *              │ +s1  +c1 │ |⟂⟩   │ +s2*c2   +(c2)² │ |⟂⟩
 *              ╰          ╯       ╰                 ╯
 */
arma::cx_mat22 OneSiteOrbitalMatrices::get_P_x_in_ge_basis(double orbital_theta) {
    const double s1 = std::sin(orbital_theta), c1 = std::cos(orbital_theta);
    const arma::mat22 P_x_in_ge_re_builder = {
        {-c1, +s1},
        {+s1, +c1}
    };
    const arma::mat22 P_x_in_ge_re = 0.5 * (arma::eye(2, 2) + P_x_in_ge_re_builder);
    const arma::cx_mat22 P_x_in_ge{P_x_in_ge_re, arma::zeros(2, 2)};
#ifndef NDEBUG
    const double s2 = std::sin(orbital_theta / 2.0), c2 = std::cos(orbital_theta / 2.0);
    const arma::mat22 P_x_in_ge_re_alternative{
        {+s2 * s2, +s2 * c2},
        {+s2 * c2, +c2 * c2}
    };
    assert(arma::norm(P_x_in_ge_re_alternative - P_x_in_ge_re) < 1e-4);
    const auto P_x_in_zx = OneSiteOrbitalMatrices::get_P_x_in_zx_basis();
    const auto beta = OneSiteOrbitalMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(P_x_in_ge - beta.t() * P_x_in_zx * beta) < 1e-4);
#endif
    return P_x_in_ge;
}


/*
 *                 |∥⟩  |⟂⟩               |∥⟩         |⟂⟩
 *              ╭          ╮           ╭                         ╮
 * P⁺ = ½ + ½ * │ +s1  +c1 │ |∥⟩ = ½ * │ 1+2*s2*c2   (c2)²-(s2)² │ |∥⟩
 *              │ +c1  -s1 │ |⟂⟩       │ (c2)²-(s2)² 1-2*s2*c2   │ |⟂⟩
 *              ╰          ╯           ╰                         ╯
 */
arma::cx_mat22 OneSiteOrbitalMatrices::get_P_plus_in_ge_basis(double orbital_theta) {
    const double s1 = std::sin(orbital_theta), c1 = std::cos(orbital_theta);
    const arma::mat22 P_plus_in_ge_re_builder = {
        {+s1, +c1},
        {+c1, -s1}
    };
    const arma::mat22 P_plus_in_ge_re = 0.5 * (arma::eye(2, 2) + P_plus_in_ge_re_builder);
    const arma::cx_mat22 P_plus_in_ge{P_plus_in_ge_re, arma::zeros(2, 2)};
#ifndef NDEBUG
    const double s2 = std::sin(orbital_theta / 2.0), c2 = std::cos(orbital_theta / 2.0);
    const arma::mat22 P_plus_in_ge_re_alternative{
        {0.5 + s2 * c2, 0.5 * (c2 * c2 - s2 * s2)},
        {0.5 * (c2 * c2 - s2 * s2), 0.5 - s2 * c2}
    };
    assert(arma::norm(P_plus_in_ge_re_alternative - P_plus_in_ge_re) < 1e-4);
    const auto P_plus_in_zx = OneSiteOrbitalMatrices::get_P_plus_in_zx_basis();
    const auto beta = OneSiteOrbitalMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(P_plus_in_ge - beta.t() * P_plus_in_zx * beta) < 1e-4);
#endif
    return P_plus_in_ge;
}


/*
 *                |∥⟩  |⟂⟩               |∥⟩         |⟂⟩
 *              ╭          ╮           ╭                         ╮
 * P⁻ = ½ + ½ * │ -s1  -c1 │ |∥⟩ = ½ * │ 1-2*s2*c2   (s2)²-(c2)² │ |∥⟩
 *              │ -c1  +s1 │ |⟂⟩       │ (s2)²-(c2)² 1+2*s2*c2   │ |⟂⟩
 *              ╰          ╯           ╰                         ╯
 */
arma::cx_mat22 OneSiteOrbitalMatrices::get_P_minus_in_ge_basis(double orbital_theta) {
    const double s1 = std::sin(orbital_theta), c1 = std::cos(orbital_theta);
    const arma::mat22 P_minus_in_ge_re_builder = {
        {-s1, -c1},
        {-c1, +s1}
    };
    const arma::mat22 P_minus_in_ge_re = 0.5 * (arma::eye(2, 2) + P_minus_in_ge_re_builder);
    const arma::cx_mat22 P_minus_in_ge{P_minus_in_ge_re, arma::zeros(2, 2)};
#ifndef NDEBUG
    const double s2 = std::sin(orbital_theta / 2.0), c2 = std::cos(orbital_theta / 2.0);
    const arma::mat22 P_minus_in_ge_re_alternative{
        {0.5 - s2 * c2, 0.5 * (s2 * s2 - c2 * c2)},
        {0.5 * (s2 * s2 - c2 * c2), 0.5 + s2 * c2}
    };
    assert(arma::norm(P_minus_in_ge_re_alternative - P_minus_in_ge_re) < 1e-4);
    const auto P_minus_in_zx = OneSiteOrbitalMatrices::get_P_minus_in_zx_basis();
    const auto beta = OneSiteOrbitalMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(P_minus_in_ge - beta.t() * P_minus_in_zx * beta) < 1e-4);
#endif
    return P_minus_in_ge;
}


arma::cx_mat22 OneSiteOrbitalMatrices::get_tau_z_in_ge_basis(double orbital_theta) {
    const arma::cx_mat22 tau_z_in_ge = OneSiteOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta) - OneSiteOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta);
#ifndef NDEBUG
    const auto tau_z_in_zx = OneSiteOrbitalMatrices::get_tau_z_in_zx_basis();
    const auto beta = OneSiteOrbitalMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(tau_z_in_ge - beta.t() * tau_z_in_zx * beta) < 1e-4);
#endif
    return OneSiteOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta) - OneSiteOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta);
}


arma::cx_mat22 OneSiteOrbitalMatrices::get_tau_x_in_ge_basis(double orbital_theta) {
    const arma::cx_mat22 tau_x_in_ge = OneSiteOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta) - OneSiteOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta);
#ifndef NDEBUG
    const auto tau_x_in_zx = OneSiteOrbitalMatrices::get_tau_x_in_zx_basis();
    const auto beta = OneSiteOrbitalMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(tau_x_in_ge - beta.t() * tau_x_in_zx * beta) < 1e-4);
#endif
    return tau_x_in_ge;
}


arma::cx_mat22 OneSiteOrbitalMatrices::get_tau_plus_in_ge_basis(double orbital_theta) {
    const arma::cx_mat22 tau_plus_in_ge = OneSiteOrbitalMatrices::get_P_plus_in_ge_basis(orbital_theta) - OneSiteOrbitalMatrices::get_P_minus_in_ge_basis(orbital_theta);
#ifndef NDEBUG
    const auto tau_plus_in_zx = OneSiteOrbitalMatrices::get_tau_plus_in_zx_basis();
    const auto beta = OneSiteOrbitalMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(tau_plus_in_ge - beta.t() * tau_plus_in_zx * beta) < 1e-4);
#endif
    return tau_plus_in_ge;
}


arma::cx_mat22 OneSiteOrbitalMatrices::get_tau_minus_in_ge_basis(double orbital_theta) {
    const arma::cx_mat22 tau_minus_in_ge = OneSiteOrbitalMatrices::get_P_minus_in_ge_basis(orbital_theta) - OneSiteOrbitalMatrices::get_P_plus_in_ge_basis(orbital_theta);
#ifndef NDEBUG
    const auto tau_minus_in_zx = OneSiteOrbitalMatrices::get_tau_minus_in_zx_basis();
    const auto beta = OneSiteOrbitalMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(tau_minus_in_ge - beta.t() * tau_minus_in_zx * beta) < 1e-4);
#endif
    return tau_minus_in_ge;
}

} // end of namespace monostar_hamiltonians

// #######################################################################
// ## OneSiteOrbitalNamedMatrices                                       ##
// #######################################################################

namespace monostar_hamiltonians {

utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_P_z_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_P_z_in_zx_basis(), "Pᶻ");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_P_x_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_P_x_in_zx_basis(), "Pˣ");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_P_plus_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_P_plus_in_zx_basis(), "P⁺");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_P_minus_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_P_minus_in_zx_basis(), "P⁻");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_tau_z_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_tau_z_in_zx_basis(), "τᶻ");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_tau_x_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_tau_x_in_zx_basis(), "τˣ");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_tau_plus_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_tau_plus_in_zx_basis(), "τ⁺");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_tau_minus_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_tau_minus_in_zx_basis(), "τ⁻");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_P_z_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta), "Pᶻ");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_P_x_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta), "Pˣ");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_P_plus_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_P_plus_in_ge_basis(orbital_theta), "P⁺");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_P_minus_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_P_minus_in_ge_basis(orbital_theta), "P⁻");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_tau_z_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_tau_z_in_ge_basis(orbital_theta), "τᶻ");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_tau_x_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_tau_x_in_ge_basis(orbital_theta), "τˣ");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_tau_plus_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_tau_plus_in_ge_basis(orbital_theta), "τ⁺");
}


utility::Named<arma::cx_mat22> OneSiteOrbitalNamedMatrices::get_tau_minus_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OneSiteOrbitalMatrices::get_tau_minus_in_ge_basis(orbital_theta), "τ⁻");
}


std::vector<utility::Named<arma::cx_mat22>> OneSiteOrbitalNamedMatrices::matrices_for_average_calculations(double orbital_theta) {
    return {
        OneSiteOrbitalNamedMatrices::get_tau_z_in_ge_basis(orbital_theta),
                OneSiteOrbitalNamedMatrices::get_tau_minus_in_ge_basis(orbital_theta),
    };
}

} // end of namespace monostar_hamiltonians

// #######################################################################
// ## TwoSitesOrbitalMatrices                                           ##
// #######################################################################

namespace monostar_hamiltonians {

arma::cx_mat44 TwoSitesOrbitalMatrices::get_P_zz_in_zx_basis() {
    const auto& P_z_in_zx_basis = OneSiteOrbitalMatrices::get_P_z_in_zx_basis();
    return arma::kron(P_z_in_zx_basis, P_z_in_zx_basis);
}


arma::cx_mat44 TwoSitesOrbitalMatrices::get_P_zx_sum_P_xz_in_zx_basis() {
    const auto& P_z_in_zx_basis = OneSiteOrbitalMatrices::get_P_z_in_zx_basis();
    const auto& P_x_in_zx_basis = OneSiteOrbitalMatrices::get_P_x_in_zx_basis();
    return arma::kron(P_z_in_zx_basis, P_x_in_zx_basis) + arma::kron(P_x_in_zx_basis, P_z_in_zx_basis);
}


arma::cx_mat44 TwoSitesOrbitalMatrices::get_P_xx_in_zx_basis() {
    const auto& P_x_in_zx_basis = OneSiteOrbitalMatrices::get_P_x_in_zx_basis();
    return arma::kron(P_x_in_zx_basis, P_x_in_zx_basis);
}


arma::cx_mat44 TwoSitesOrbitalMatrices::get_P_zz_in_ge_basis(double orbital_theta) {
    const auto& P_z_in_ge_basis = OneSiteOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta);
    return arma::kron(P_z_in_ge_basis, P_z_in_ge_basis);
}


arma::cx_mat44 TwoSitesOrbitalMatrices::get_P_zx_sum_P_xz_in_ge_basis(double orbital_theta) {
    const auto& P_z_in_ge_basis = OneSiteOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta);
    const auto& P_x_in_ge_basis = OneSiteOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta);
    return arma::kron(P_z_in_ge_basis, P_x_in_ge_basis) + arma::kron(P_x_in_ge_basis, P_z_in_ge_basis);
}


arma::cx_mat44 TwoSitesOrbitalMatrices::get_P_xx_in_ge_basis(double orbital_theta) {
    const auto& P_x_in_ge_basis = OneSiteOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta);
    return arma::kron(P_x_in_ge_basis, P_x_in_ge_basis);
}

} // end of namespace monostar_hamiltonians

// #######################################################################
// ## TwoSitesOrbitalNamedMatrices                                      ##
// #######################################################################

namespace monostar_hamiltonians {

utility::Named<arma::cx_mat44> TwoSitesOrbitalNamedMatrices::get_P_zz_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(TwoSitesOrbitalMatrices::get_P_zz_in_zx_basis(), "Pᶻᶻ");
}


utility::Named<arma::cx_mat44> TwoSitesOrbitalNamedMatrices::get_P_zx_sum_P_xz_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(TwoSitesOrbitalMatrices::get_P_zx_sum_P_xz_in_zx_basis(), "(Pᶻˣ+Pˣᶻ)");
}


utility::Named<arma::cx_mat44> TwoSitesOrbitalNamedMatrices::get_P_xx_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(TwoSitesOrbitalMatrices::get_P_xx_in_zx_basis(), "Pˣˣ");
}


utility::Named<arma::cx_mat44> TwoSitesOrbitalNamedMatrices::get_P_zz_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(TwoSitesOrbitalMatrices::get_P_zz_in_ge_basis(orbital_theta), "Pᶻᶻ");
}


utility::Named<arma::cx_mat44> TwoSitesOrbitalNamedMatrices::get_P_zx_sum_P_xz_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(TwoSitesOrbitalMatrices::get_P_zx_sum_P_xz_in_ge_basis(orbital_theta), "(Pᶻˣ+Pˣᶻ)");
}


utility::Named<arma::cx_mat44> TwoSitesOrbitalNamedMatrices::get_P_xx_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(TwoSitesOrbitalMatrices::get_P_xx_in_ge_basis(orbital_theta), "Pˣˣ");
}


std::vector<utility::Named<arma::cx_mat44>> TwoSitesOrbitalNamedMatrices::matrices_for_average_calculations(double orbital_theta) {
    return {
        TwoSitesOrbitalNamedMatrices::get_P_zz_in_ge_basis(orbital_theta),
                TwoSitesOrbitalNamedMatrices::get_P_zx_sum_P_xz_in_ge_basis(orbital_theta),
                TwoSitesOrbitalNamedMatrices::get_P_xx_in_ge_basis(orbital_theta)
    };
}

} // end of namespace monostar_hamiltonians
