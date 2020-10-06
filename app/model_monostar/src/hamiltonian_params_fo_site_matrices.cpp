#include<model_monostar/hamiltonian_params_fo_site_matrices.hpp>

#include<cmath>
#include<cassert>


// #######################################################################
// ## OrbitalSiteMatrices                                               ##
// #######################################################################

arma::cx_mat22 OrbitalSiteMatrices::get_P_z_in_zx_basis() {
    const arma::mat22 P_z_in_zx_basis_re = {
        {1.0, 0.0},
        {0.0, 0.0}};
    arma::cx_mat P_z_in_zx_basis{P_z_in_zx_basis_re, arma::zeros(2, 2)};
    return P_z_in_zx_basis;
}


arma::cx_mat22 OrbitalSiteMatrices::get_P_x_in_zx_basis() {
    const arma::mat22 P_x_in_zx_basis_re = {
        {0.0, 0.0},
        {0.0, 1.0}};
    arma::cx_mat P_x_in_zx_basis{P_x_in_zx_basis_re, arma::zeros(2, 2)};
    return P_x_in_zx_basis;
}


arma::cx_mat22 OrbitalSiteMatrices::get_P_plus_in_zx_basis() {
    const arma::mat22 P_plus_in_zx_basis_re = {
        {+0.5, +0.5},
        {+0.5, +0.5}};
    arma::cx_mat P_plus_in_zx_basis{P_plus_in_zx_basis_re, arma::zeros(2, 2)};
    return P_plus_in_zx_basis;
}


arma::cx_mat22 OrbitalSiteMatrices::get_P_minus_in_zx_basis() {
    const arma::mat22 P_minus_in_zx_basis_re = {
        {+0.5, -0.5},
        {-0.5, +0.5}};
    arma::cx_mat P_minus_in_zx_basis{P_minus_in_zx_basis_re, arma::zeros(2, 2)};
    return P_minus_in_zx_basis;
}


arma::cx_mat22 OrbitalSiteMatrices::get_tau_z_in_zx_basis() {
    return OrbitalSiteMatrices::get_P_z_in_zx_basis() - OrbitalSiteMatrices::get_P_x_in_zx_basis();
}


arma::cx_mat22 OrbitalSiteMatrices::get_tau_x_in_zx_basis() {
    return OrbitalSiteMatrices::get_P_x_in_zx_basis() - OrbitalSiteMatrices::get_P_z_in_zx_basis();
}


arma::cx_mat22 OrbitalSiteMatrices::get_tau_plus_in_zx_basis() {
    return OrbitalSiteMatrices::get_P_plus_in_zx_basis() - OrbitalSiteMatrices::get_P_minus_in_zx_basis();
}


arma::cx_mat22 OrbitalSiteMatrices::get_tau_minus_in_zx_basis() {
    return OrbitalSiteMatrices::get_P_minus_in_zx_basis() - OrbitalSiteMatrices::get_P_plus_in_zx_basis();
}


/*
 * ╭         ╮     ╭         ╮              ╭                        ╮
 * │ |∥⟩ |⟂⟩ │  =  │ |x⟩ |z⟩ │ β, where β = │  +cos(θ/2)   -sin(θ/2) │
 * ╰         ╯     ╰         ╯              │  +sin(θ/2)   +cos(θ/2) │
 *                                          ╰                        ╯
 */
arma::cx_mat22 OrbitalSiteMatrices::get_beta_from_zx_to_ge(double orbital_theta) {
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
arma::cx_mat22 OrbitalSiteMatrices::get_P_z_in_ge_basis(double orbital_theta) {
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
    const auto P_z_in_zx = OrbitalSiteMatrices::get_P_z_in_zx_basis();
    const auto beta = OrbitalSiteMatrices::get_beta_from_zx_to_ge(orbital_theta);
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
arma::cx_mat22 OrbitalSiteMatrices::get_P_x_in_ge_basis(double orbital_theta) {
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
    const auto P_x_in_zx = OrbitalSiteMatrices::get_P_x_in_zx_basis();
    const auto beta = OrbitalSiteMatrices::get_beta_from_zx_to_ge(orbital_theta);
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
arma::cx_mat22 OrbitalSiteMatrices::get_P_plus_in_ge_basis(double orbital_theta) {
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
    const auto P_plus_in_zx = OrbitalSiteMatrices::get_P_plus_in_zx_basis();
    const auto beta = OrbitalSiteMatrices::get_beta_from_zx_to_ge(orbital_theta);
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
arma::cx_mat22 OrbitalSiteMatrices::get_P_minus_in_ge_basis(double orbital_theta) {
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
    const auto P_minus_in_zx = OrbitalSiteMatrices::get_P_minus_in_zx_basis();
    const auto beta = OrbitalSiteMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(P_minus_in_ge - beta.t() * P_minus_in_zx * beta) < 1e-4);
#endif
    return P_minus_in_ge;
}


arma::cx_mat22 OrbitalSiteMatrices::get_tau_z_in_ge_basis(double orbital_theta) {
    const arma::cx_mat22 tau_z_in_ge = OrbitalSiteMatrices::get_P_z_in_ge_basis(orbital_theta) - OrbitalSiteMatrices::get_P_x_in_ge_basis(orbital_theta);
#ifndef NDEBUG
    const auto tau_z_in_zx = OrbitalSiteMatrices::get_tau_z_in_zx_basis();
    const auto beta = OrbitalSiteMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(tau_z_in_ge - beta.t() * tau_z_in_zx * beta) < 1e-4);
#endif
    return OrbitalSiteMatrices::get_P_z_in_ge_basis(orbital_theta) - OrbitalSiteMatrices::get_P_x_in_ge_basis(orbital_theta);
}


arma::cx_mat22 OrbitalSiteMatrices::get_tau_x_in_ge_basis(double orbital_theta) {
    const arma::cx_mat22 tau_x_in_ge = OrbitalSiteMatrices::get_P_x_in_ge_basis(orbital_theta) - OrbitalSiteMatrices::get_P_z_in_ge_basis(orbital_theta);
#ifndef NDEBUG
    const auto tau_x_in_zx = OrbitalSiteMatrices::get_tau_x_in_zx_basis();
    const auto beta = OrbitalSiteMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(tau_x_in_ge - beta.t() * tau_x_in_zx * beta) < 1e-4);
#endif
    return tau_x_in_ge;
}


arma::cx_mat22 OrbitalSiteMatrices::get_tau_plus_in_ge_basis(double orbital_theta) {
    const arma::cx_mat22 tau_plus_in_ge = OrbitalSiteMatrices::get_P_plus_in_ge_basis(orbital_theta) - OrbitalSiteMatrices::get_P_minus_in_ge_basis(orbital_theta);
#ifndef NDEBUG
    const auto tau_plus_in_zx = OrbitalSiteMatrices::get_tau_plus_in_zx_basis();
    const auto beta = OrbitalSiteMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(tau_plus_in_ge - beta.t() * tau_plus_in_zx * beta) < 1e-4);
#endif
    return tau_plus_in_ge;
}


arma::cx_mat22 OrbitalSiteMatrices::get_tau_minus_in_ge_basis(double orbital_theta) {
    const arma::cx_mat22 tau_minus_in_ge = OrbitalSiteMatrices::get_P_minus_in_ge_basis(orbital_theta) - OrbitalSiteMatrices::get_P_plus_in_ge_basis(orbital_theta);
#ifndef NDEBUG
    const auto tau_minus_in_zx = OrbitalSiteMatrices::get_tau_minus_in_zx_basis();
    const auto beta = OrbitalSiteMatrices::get_beta_from_zx_to_ge(orbital_theta);
    assert(arma::norm(tau_minus_in_ge - beta.t() * tau_minus_in_zx * beta) < 1e-4);
#endif
    return tau_minus_in_ge;
}


// #######################################################################
// ## OrbitalSiteNamedMatrices                                          ##
// #######################################################################


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_P_z_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_P_z_in_zx_basis(), "Pᶻ_in_zx");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_P_x_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_P_x_in_zx_basis(), "Pˣ_in_zx");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_P_plus_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_P_plus_in_zx_basis(), "P⁺ [in zx]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_P_minus_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_P_minus_in_zx_basis(), "P⁻ [in zx]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_tau_z_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_tau_z_in_zx_basis(), "τᶻ [in zx]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_tau_x_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_tau_x_in_zx_basis(), "τˣ [in zx]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_tau_plus_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_tau_plus_in_zx_basis(), "τ⁺ [in zx]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_tau_minus_in_zx_basis() {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_tau_minus_in_zx_basis(), "τ⁻ [in zx]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_P_z_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_P_z_in_ge_basis(orbital_theta), "Pᶻ [in ge]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_P_x_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_P_x_in_ge_basis(orbital_theta), "Pˣ[in ge]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_P_plus_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_P_plus_in_ge_basis(orbital_theta), "P⁺ [in_ge]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_P_minus_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_P_minus_in_ge_basis(orbital_theta), "P⁻ [in ge]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_tau_z_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_tau_z_in_ge_basis(orbital_theta), "τᶻ [in ge]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_tau_x_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_tau_x_in_ge_basis(orbital_theta), "τˣ [in ge]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_tau_plus_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_tau_plus_in_ge_basis(orbital_theta), "τ⁺ [in ge]");
}


utility::Named<arma::cx_mat22> OrbitalSiteNamedMatrices::get_tau_minus_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat22>(OrbitalSiteMatrices::get_tau_minus_in_ge_basis(orbital_theta), "τ⁻ [in ge]");
}

// #######################################################################
// ## OrbitalSiteNamedMatricesAverageCalculations                       ##
// #######################################################################

std::vector<utility::Named<arma::cx_mat22>> orbital_site_matrices_for_average_calculations(double orbital_theta) {
    return {
        OrbitalSiteNamedMatrices::get_tau_z_in_ge_basis(orbital_theta),
                OrbitalSiteNamedMatrices::get_tau_minus_in_ge_basis(orbital_theta),
    };
}

// #######################################################################
// ## OrbitalTwoSiteMatrices                                            ##
// #######################################################################

arma::cx_mat44 OrbitalTwoSiteMatrices::get_P_zz_in_zx_basis() {
    const auto& P_z_in_zx_basis = OrbitalSiteMatrices::get_P_z_in_zx_basis();
    return arma::kron(P_z_in_zx_basis, P_z_in_zx_basis);
}


arma::cx_mat44 OrbitalTwoSiteMatrices::get_P_zx_sum_P_xz_in_zx_basis() {
    const auto& P_z_in_zx_basis = OrbitalSiteMatrices::get_P_z_in_zx_basis();
    const auto& P_x_in_zx_basis = OrbitalSiteMatrices::get_P_x_in_zx_basis();
    return arma::kron(P_z_in_zx_basis, P_x_in_zx_basis) + arma::kron(P_x_in_zx_basis, P_z_in_zx_basis);
}


arma::cx_mat44 OrbitalTwoSiteMatrices::get_P_xx_in_zx_basis() {
    const auto& P_x_in_zx_basis = OrbitalSiteMatrices::get_P_x_in_zx_basis();
    return arma::kron(P_x_in_zx_basis, P_x_in_zx_basis);
}


arma::cx_mat44 OrbitalTwoSiteMatrices::get_P_zz_in_ge_basis(double orbital_theta) {
    const auto& P_z_in_ge_basis = OrbitalSiteMatrices::get_P_z_in_ge_basis(orbital_theta);
    return arma::kron(P_z_in_ge_basis, P_z_in_ge_basis);
}


arma::cx_mat44 OrbitalTwoSiteMatrices::get_P_zx_sum_P_xz_in_ge_basis(double orbital_theta) {
    const auto& P_z_in_ge_basis = OrbitalSiteMatrices::get_P_z_in_ge_basis(orbital_theta);
    const auto& P_x_in_ge_basis = OrbitalSiteMatrices::get_P_x_in_ge_basis(orbital_theta);
    return arma::kron(P_z_in_ge_basis, P_x_in_ge_basis) + arma::kron(P_x_in_ge_basis, P_z_in_ge_basis);
}


arma::cx_mat44 OrbitalTwoSiteMatrices::get_P_xx_in_ge_basis(double orbital_theta) {
    const auto& P_x_in_ge_basis = OrbitalSiteMatrices::get_P_x_in_ge_basis(orbital_theta);
    return arma::kron(P_x_in_ge_basis, P_x_in_ge_basis);
}

// #######################################################################
// ## OrbitalTwoSiteNamedMatrices                                       ##
// #######################################################################

utility::Named<arma::cx_mat44> OrbitalTwoSiteNamedMatrices::get_P_zz_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(OrbitalTwoSiteMatrices::get_P_zz_in_zx_basis(), "Pᶻᶻ [in zx]");
}


utility::Named<arma::cx_mat44> OrbitalTwoSiteNamedMatrices::get_P_zx_sum_P_xz_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(OrbitalTwoSiteMatrices::get_P_zx_sum_P_xz_in_zx_basis(), "(Pᶻˣ+Pˣᶻ) [in zx]");
}


utility::Named<arma::cx_mat44> OrbitalTwoSiteNamedMatrices::get_P_xx_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(OrbitalTwoSiteMatrices::get_P_xx_in_zx_basis(), "Pˣˣ [in_zx]");
}


utility::Named<arma::cx_mat44> OrbitalTwoSiteNamedMatrices::get_P_zz_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(OrbitalTwoSiteMatrices::get_P_zz_in_ge_basis(orbital_theta), "Pᶻᶻ [in ge]");
}


utility::Named<arma::cx_mat44> OrbitalTwoSiteNamedMatrices::get_P_zx_sum_P_xz_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(OrbitalTwoSiteMatrices::get_P_zx_sum_P_xz_in_ge_basis(orbital_theta), "(Pᶻˣ+Pˣᶻ) [in ge]");
}


utility::Named<arma::cx_mat44> OrbitalTwoSiteNamedMatrices::get_P_xx_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(OrbitalTwoSiteMatrices::get_P_xx_in_ge_basis(orbital_theta), "Pˣˣ [in ge]");
}

// #######################################################################
// ## OrbitalTwoSiteNamedMatricesForAverageCalculations                 ##
// #######################################################################

std::vector<utility::Named<arma::cx_mat44>> orbital_two_site_matrices_for_average_calculations(double orbital_theta) {
    return {
        OrbitalTwoSiteNamedMatrices::get_P_zz_in_ge_basis(orbital_theta),
                OrbitalTwoSiteNamedMatrices::get_P_zx_sum_P_xz_in_ge_basis(orbital_theta),
                OrbitalTwoSiteNamedMatrices::get_P_xx_in_ge_basis(orbital_theta)
    };
}
