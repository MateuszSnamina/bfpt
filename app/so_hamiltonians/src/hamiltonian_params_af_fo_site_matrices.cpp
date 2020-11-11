#include <so_hamiltonians/hamiltonian_params_af_fo_site_matrices.hpp>

#include <monostar_hamiltonians/hamiltonian_params_fo_site_matrices.hpp>
#include <monostar_hamiltonians/hamiltonian_params_af_fm_site_matrices.hpp>

#include <cmath>
#include <cassert>

// #######################################################################
// ## OneSiteSpinOrbitalMatrices                                        ##
// #######################################################################

namespace so_hamiltonians {

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_S_z_af_latticeA() {
    return arma::kron(
        monostar_hamiltonians::OneSiteSpinMatrices::get_S_z_in_ge_af_latticeA(),
        arma::eye(2, 2));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_S_z_af_latticeB() {
    return arma::kron(
        monostar_hamiltonians::OneSiteSpinMatrices::get_S_z_in_ge_af_latticeB(),
        arma::eye(2, 2));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_S_plus_af_latticeA() {
    return arma::kron(
        monostar_hamiltonians::OneSiteSpinMatrices::get_S_plus_in_ge_af_latticeA(),
        arma::eye(2, 2));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_S_plus_af_latticeB() {
    return arma::kron(
        monostar_hamiltonians::OneSiteSpinMatrices::get_S_plus_in_ge_af_latticeB(),
        arma::eye(2, 2));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_S_minus_af_latticeA() {
    return arma::kron(
        monostar_hamiltonians::OneSiteSpinMatrices::get_S_minus_in_ge_af_latticeA(),
        arma::eye(2, 2));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_S_minus_af_latticeB() {
    return arma::kron(
        monostar_hamiltonians::OneSiteSpinMatrices::get_S_minus_in_ge_af_latticeB(),
        arma::eye(2, 2));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_P_z_in_zx_basis() {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_z_in_zx_basis());
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_P_x_in_zx_basis() {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_x_in_zx_basis());
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_P_plus_in_zx_basis() {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_plus_in_zx_basis());
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_P_minus_in_zx_basis() {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_minus_in_zx_basis());
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_tau_z_in_zx_basis() {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_z_in_zx_basis());
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_tau_x_in_zx_basis() {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_x_in_zx_basis());
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_tau_plus_in_zx_basis() {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_plus_in_zx_basis());
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_tau_minus_in_zx_basis() {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_minus_in_zx_basis());
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_beta_from_zx_to_ge(double orbital_theta) {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_beta_from_zx_to_ge(orbital_theta));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_P_z_in_ge_basis(double orbital_theta) {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_P_x_in_ge_basis(double orbital_theta) {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_P_plus_in_ge_basis(double orbital_theta) {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_plus_in_ge_basis(orbital_theta));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_P_minus_in_ge_basis(double orbital_theta) {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_minus_in_ge_basis(orbital_theta));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_tau_z_in_ge_basis(double orbital_theta) {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_z_in_ge_basis(orbital_theta));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_tau_x_in_ge_basis(double orbital_theta) {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_x_in_ge_basis(orbital_theta));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_tau_plus_in_ge_basis(double orbital_theta) {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_plus_in_ge_basis(orbital_theta));
}

arma::cx_mat44 OneSiteSpinOrbitalMatrices::get_tau_minus_in_ge_basis(double orbital_theta) {
    return arma::kron(
        arma::eye(2, 2),
        monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_minus_in_ge_basis(orbital_theta));
}

}  // end of namespace so_hamiltonians

// #######################################################################
// ## OneSiteSpinOrbitalNamedMatrices                                   ##
// #######################################################################

namespace so_hamiltonians {

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_S_z_af_latticeA() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_S_z_af_latticeA(), "Sᶻ@A");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_S_z_af_latticeB() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_S_z_af_latticeB(), "Sᶻ@B");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_S_plus_af_latticeA() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_S_plus_af_latticeA(), "S⁺@A");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_S_plus_af_latticeB() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_S_plus_af_latticeB(), "S⁺@B");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_S_minus_af_latticeA() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_S_minus_af_latticeA(), "S⁻@A");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_S_minus_af_latticeB() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_S_minus_af_latticeB(), "S⁻@B");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_P_z_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_P_z_in_zx_basis(), "Pᶻ");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_P_x_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_P_x_in_zx_basis(), "Pˣ");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_P_plus_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_P_plus_in_zx_basis(), "P⁺");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_P_minus_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_P_minus_in_zx_basis(), "P⁻");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_tau_z_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_tau_z_in_zx_basis(), "τᶻ");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_tau_x_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_tau_x_in_zx_basis(), "τˣ");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_tau_plus_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_tau_plus_in_zx_basis(), "τ⁺");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_tau_minus_in_zx_basis() {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_tau_minus_in_zx_basis(), "τ⁻");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_P_z_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta), "Pᶻ");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_P_x_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta), "Pˣ");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_P_plus_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_P_plus_in_ge_basis(orbital_theta), "P⁺");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_P_minus_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_P_minus_in_ge_basis(orbital_theta), "P⁻");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_tau_z_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_tau_z_in_ge_basis(orbital_theta), "τᶻ");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_tau_x_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_tau_x_in_ge_basis(orbital_theta), "τˣ");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_tau_plus_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_tau_plus_in_ge_basis(orbital_theta), "τ⁺");
}

utility::Named<arma::cx_mat44> OneSiteSpinOrbitalNamedMatrices::get_tau_minus_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::cx_mat44>(OneSiteSpinOrbitalMatrices::get_tau_minus_in_ge_basis(orbital_theta), "τ⁻");
}

std::vector<utility::Named<arma::cx_mat44>> OneSiteSpinOrbitalNamedMatrices::matrices_for_average_calculations(double orbital_theta) {
    return {
        OneSiteSpinOrbitalNamedMatrices::get_tau_z_in_ge_basis(orbital_theta),
        OneSiteSpinOrbitalNamedMatrices::get_tau_minus_in_ge_basis(orbital_theta),
    };
}

}  // end of namespace so_hamiltonians

// #######################################################################
// ## TwoSitesSpinOrbitalMatrices                                       ##
// #######################################################################

namespace so_hamiltonians {

arma::extension::cx_mat1616 TwoSitesSpinOrbitalMatrices::get_S_S_ondiag_af() {
    const arma::extension::cx_mat1616 S_S_ondiag_af_latticeAB =
        arma::kron(OneSiteSpinOrbitalMatrices::get_S_z_af_latticeA(), OneSiteSpinOrbitalMatrices::get_S_z_af_latticeB());
#ifndef NDEBUG
    const arma::extension::cx_mat1616 S_S_ondiag_af_latticeBA =
        arma::kron(OneSiteSpinOrbitalMatrices::get_S_z_af_latticeB(), OneSiteSpinOrbitalMatrices::get_S_z_af_latticeA());
    assert(arma::norm(S_S_ondiag_af_latticeAB - S_S_ondiag_af_latticeBA) < 1e-8);
#endif
    return S_S_ondiag_af_latticeAB;
}

arma::extension::cx_mat1616 TwoSitesSpinOrbitalMatrices::get_S_S_offdiag_af() {
    const arma::extension::cx_mat1616 S_S_offdiagaf_latticeAB = 0.5 * (arma::kron(OneSiteSpinOrbitalMatrices::get_S_plus_af_latticeA(), OneSiteSpinOrbitalMatrices::get_S_minus_af_latticeB()) +
                                                                       arma::kron(OneSiteSpinOrbitalMatrices::get_S_minus_af_latticeA(), OneSiteSpinOrbitalMatrices::get_S_plus_af_latticeB()));
#ifndef NDEBUG
    const arma::extension::cx_mat1616 S_S_offdiag_af_latticeBA = 0.5 * (arma::kron(OneSiteSpinOrbitalMatrices::get_S_plus_af_latticeB(), OneSiteSpinOrbitalMatrices::get_S_minus_af_latticeA()) +
                                                                        arma::kron(OneSiteSpinOrbitalMatrices::get_S_minus_af_latticeB(), OneSiteSpinOrbitalMatrices::get_S_plus_af_latticeA()));
    assert(arma::norm(S_S_offdiagaf_latticeAB - S_S_offdiag_af_latticeBA) < 1e-8);
#endif
    return S_S_offdiagaf_latticeAB;
}

arma::extension::cx_mat1616 TwoSitesSpinOrbitalMatrices::get_S_S_af() {
    return get_S_S_ondiag_af() + get_S_S_offdiag_af();
}

arma::extension::cx_mat1616 TwoSitesSpinOrbitalMatrices::get_P_zz_in_zx_basis() {
    const auto& P_z_in_zx_basis = OneSiteSpinOrbitalMatrices::get_P_z_in_zx_basis();
    return arma::kron(P_z_in_zx_basis, P_z_in_zx_basis);
}

arma::extension::cx_mat1616 TwoSitesSpinOrbitalMatrices::get_P_zx_sum_P_xz_in_zx_basis() {
    const auto& P_z_in_zx_basis = OneSiteSpinOrbitalMatrices::get_P_z_in_zx_basis();
    const auto& P_x_in_zx_basis = OneSiteSpinOrbitalMatrices::get_P_x_in_zx_basis();
    return arma::kron(P_z_in_zx_basis, P_x_in_zx_basis) + arma::kron(P_x_in_zx_basis, P_z_in_zx_basis);
}

arma::extension::cx_mat1616 TwoSitesSpinOrbitalMatrices::get_P_xx_in_zx_basis() {
    const auto& P_x_in_zx_basis = OneSiteSpinOrbitalMatrices::get_P_x_in_zx_basis();
    return arma::kron(P_x_in_zx_basis, P_x_in_zx_basis);
}

arma::extension::cx_mat1616 TwoSitesSpinOrbitalMatrices::get_P_zz_in_ge_basis(double orbital_theta) {
    const auto& P_z_in_ge_basis = OneSiteSpinOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta);
    return arma::kron(P_z_in_ge_basis, P_z_in_ge_basis);
}

arma::extension::cx_mat1616 TwoSitesSpinOrbitalMatrices::get_P_zx_sum_P_xz_in_ge_basis(double orbital_theta) {
    const auto& P_z_in_ge_basis = OneSiteSpinOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta);
    const auto& P_x_in_ge_basis = OneSiteSpinOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta);
    return arma::kron(P_z_in_ge_basis, P_x_in_ge_basis) + arma::kron(P_x_in_ge_basis, P_z_in_ge_basis);
}

arma::extension::cx_mat1616 TwoSitesSpinOrbitalMatrices::get_P_xx_in_ge_basis(double orbital_theta) {
    const auto& P_x_in_ge_basis = OneSiteSpinOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta);
    return arma::kron(P_x_in_ge_basis, P_x_in_ge_basis);
}

}  // end of namespace so_hamiltonians

// #######################################################################
// ## TwoSitesSpinOrbitalNamedMatrices                                  ##
// #######################################################################

namespace so_hamiltonians {

utility::Named<arma::extension::cx_mat1616> TwoSitesSpinOrbitalNamedMatrices::get_S_S_ondiag_af() {
    return utility::Named<arma::extension::cx_mat1616>(TwoSitesSpinOrbitalMatrices::get_S_S_ondiag_af(), "SᶻSᶻ");
}

utility::Named<arma::extension::cx_mat1616> TwoSitesSpinOrbitalNamedMatrices::get_S_S_offdiag_af() {
    return utility::Named<arma::extension::cx_mat1616>(TwoSitesSpinOrbitalMatrices::get_S_S_offdiag_af(), "½(S⁺S⁻+S⁻S⁺)");
}

utility::Named<arma::extension::cx_mat1616> TwoSitesSpinOrbitalNamedMatrices::get_S_S_af() {
    return utility::Named<arma::extension::cx_mat1616>(TwoSitesSpinOrbitalMatrices::get_S_S_af(), "S⋅S");
}

utility::Named<arma::extension::cx_mat1616> TwoSitesSpinOrbitalNamedMatrices::get_P_zz_in_zx_basis() {
    return utility::Named<arma::extension::cx_mat1616>(TwoSitesSpinOrbitalMatrices::get_P_zz_in_zx_basis(), "Pᶻᶻ");
}

utility::Named<arma::extension::cx_mat1616> TwoSitesSpinOrbitalNamedMatrices::get_P_zx_sum_P_xz_in_zx_basis() {
    return utility::Named<arma::extension::cx_mat1616>(TwoSitesSpinOrbitalMatrices::get_P_zx_sum_P_xz_in_zx_basis(), "(Pᶻˣ+Pˣᶻ)");
}

utility::Named<arma::extension::cx_mat1616> TwoSitesSpinOrbitalNamedMatrices::get_P_xx_in_zx_basis() {
    return utility::Named<arma::extension::cx_mat1616>(TwoSitesSpinOrbitalMatrices::get_P_xx_in_zx_basis(), "Pˣˣ");
}

utility::Named<arma::extension::cx_mat1616> TwoSitesSpinOrbitalNamedMatrices::get_P_zz_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::extension::cx_mat1616>(TwoSitesSpinOrbitalMatrices::get_P_zz_in_ge_basis(orbital_theta), "Pᶻᶻ");
}

utility::Named<arma::extension::cx_mat1616> TwoSitesSpinOrbitalNamedMatrices::get_P_zx_sum_P_xz_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::extension::cx_mat1616>(TwoSitesSpinOrbitalMatrices::get_P_zx_sum_P_xz_in_ge_basis(orbital_theta), "Pᶻˣ+Pˣᶻ");
}

utility::Named<arma::extension::cx_mat1616> TwoSitesSpinOrbitalNamedMatrices::get_P_xx_in_ge_basis(double orbital_theta) {
    return utility::Named<arma::extension::cx_mat1616>(TwoSitesSpinOrbitalMatrices::get_P_xx_in_ge_basis(orbital_theta), "Pˣˣ");
}

std::vector<utility::Named<arma::extension::cx_mat1616>> TwoSitesSpinOrbitalNamedMatrices::matrices_for_average_calculations(double orbital_theta) {
    return {
        TwoSitesSpinOrbitalNamedMatrices::get_S_S_ondiag_af(),
        TwoSitesSpinOrbitalNamedMatrices::get_S_S_offdiag_af(),
        TwoSitesSpinOrbitalNamedMatrices::get_S_S_af(),
        TwoSitesSpinOrbitalNamedMatrices::get_P_zz_in_ge_basis(orbital_theta),
        TwoSitesSpinOrbitalNamedMatrices::get_P_zx_sum_P_xz_in_ge_basis(orbital_theta),
        TwoSitesSpinOrbitalNamedMatrices::get_P_xx_in_ge_basis(orbital_theta)};
}

}  // end of namespace so_hamiltonians
