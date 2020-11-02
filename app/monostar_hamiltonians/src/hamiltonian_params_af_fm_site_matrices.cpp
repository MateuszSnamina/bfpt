#include <monostar_hamiltonians/hamiltonian_params_af_fm_site_matrices.hpp>

#include <cassert>

// #######################################################################
// ## Helpers -- Pauli matrices                                         ##
// #######################################################################

namespace {

struct PaluliMatrices {
    static arma::cx_mat22 get_sigma_z();
    static arma::cx_mat22 get_sigma_plus();
    static arma::cx_mat22 get_sigma_minus();
};

arma::cx_mat22 PaluliMatrices::get_sigma_z() {
    arma::mat22 sigma_z_re{
        {+1, 0},
        {0, -1}};
    return arma::cx_mat22(sigma_z_re, arma::zeros(2, 2));
}

arma::cx_mat22 PaluliMatrices::get_sigma_plus() {
    arma::mat22 sigma_plus_re{
        {0, 1},
        {0, 0}};
    return arma::cx_mat22(2 * sigma_plus_re, arma::zeros(2, 2));
}

arma::cx_mat22 PaluliMatrices::get_sigma_minus() {
    arma::mat22 sigma_minus_re{
        {0, 0},
        {1, 0}};
    return arma::cx_mat22(2 * sigma_minus_re, arma::zeros(2, 2));
}

}  // end of namespace

// #######################################################################
// ## OneSiteSpinMatrices                                               ##
// #######################################################################

namespace monostar_hamiltonians {

arma::cx_mat22 OneSiteSpinMatrices::get_S_z_fm() {
    return 0.5 * PaluliMatrices::get_sigma_z();
}

arma::cx_mat22 OneSiteSpinMatrices::get_S_z_af_latticeA() {
    return +0.5 * PaluliMatrices::get_sigma_z();
}

arma::cx_mat22 OneSiteSpinMatrices::get_S_z_af_latticeB() {
    return -0.5 * PaluliMatrices::get_sigma_z();
}

arma::cx_mat22 OneSiteSpinMatrices::get_S_plus_fm() {
    return 0.5 * PaluliMatrices::get_sigma_plus();
}

arma::cx_mat22 OneSiteSpinMatrices::get_S_plus_af_latticeA() {
    return 0.5 * PaluliMatrices::get_sigma_plus();
}

arma::cx_mat22 OneSiteSpinMatrices::get_S_plus_af_latticeB() {
    return 0.5 * PaluliMatrices::get_sigma_minus();
}

arma::cx_mat22 OneSiteSpinMatrices::get_S_minus_fm() {
    return 0.5 * PaluliMatrices::get_sigma_minus();
}

arma::cx_mat22 OneSiteSpinMatrices::get_S_minus_af_latticeA() {
    return 0.5 * PaluliMatrices::get_sigma_minus();
}

arma::cx_mat22 OneSiteSpinMatrices::get_S_minus_af_latticeB() {
    return 0.5 * PaluliMatrices::get_sigma_plus();
}

}  // end of namespace monostar_hamiltonians

// #######################################################################
// ## OneSiteSpinNamedMatrices                                          ##
// #######################################################################

namespace monostar_hamiltonians {

utility::Named<arma::cx_mat22> OneSiteSpinNamedMatrices::get_S_z_fm() {
    return utility::Named<arma::cx_mat22>(OneSiteSpinMatrices::get_S_z_fm(), "Sᶻ");
}

utility::Named<arma::cx_mat22> OneSiteSpinNamedMatrices::get_S_z_af_latticeA() {
    return utility::Named<arma::cx_mat22>(OneSiteSpinMatrices::get_S_z_af_latticeA(), "Sᶻ@A");
}

utility::Named<arma::cx_mat22> OneSiteSpinNamedMatrices::get_S_z_af_latticeB() {
    return utility::Named<arma::cx_mat22>(OneSiteSpinMatrices::get_S_z_af_latticeB(), "Sᶻ@B");
}

utility::Named<arma::cx_mat22> OneSiteSpinNamedMatrices::get_S_plus_fm() {
    return utility::Named<arma::cx_mat22>(OneSiteSpinMatrices::get_S_plus_fm(), "S⁺");
}

utility::Named<arma::cx_mat22> OneSiteSpinNamedMatrices::get_S_plus_af_latticeA() {
    return utility::Named<arma::cx_mat22>(OneSiteSpinMatrices::get_S_plus_af_latticeA(), "S⁺@A");
}

utility::Named<arma::cx_mat22> OneSiteSpinNamedMatrices::get_S_plus_af_latticeB() {
    return utility::Named<arma::cx_mat22>(OneSiteSpinMatrices::get_S_plus_af_latticeB(), "S⁺@B");
}

utility::Named<arma::cx_mat22> OneSiteSpinNamedMatrices::get_S_minus_fm() {
    return utility::Named<arma::cx_mat22>(OneSiteSpinMatrices::get_S_minus_fm(), "S⁻");
}

utility::Named<arma::cx_mat22> OneSiteSpinNamedMatrices::get_S_minus_af_latticeA() {
    return utility::Named<arma::cx_mat22>(OneSiteSpinMatrices::get_S_minus_af_latticeA(), "S⁻@A");
}

utility::Named<arma::cx_mat22> OneSiteSpinNamedMatrices::get_S_minus_af_latticeB() {
    return utility::Named<arma::cx_mat22>(OneSiteSpinMatrices::get_S_minus_af_latticeB(), "S⁻@B");
}

std::vector<utility::Named<arma::cx_mat22>> OneSiteSpinNamedMatrices::site_matrices_for_average_calculations_af_fm() {
    return {};
}

}  // end of namespace monostar_hamiltonians

// #######################################################################
// ## TwoSitesSpinNamedMatrices                                         ##
// #######################################################################

namespace monostar_hamiltonians {

arma::cx_mat44 TwoSitesSpinMatrices::get_S_S_ondiag_fm() {
    return arma::kron(OneSiteSpinMatrices::get_S_z_fm(), OneSiteSpinMatrices::get_S_z_fm());
}

arma::cx_mat44 TwoSitesSpinMatrices::get_S_S_ondiag_af() {
    const arma::cx_mat44 S_S_ondiag_af_latticeAB = arma::kron(OneSiteSpinMatrices::get_S_z_af_latticeA(), OneSiteSpinMatrices::get_S_z_af_latticeB());
#ifndef NDEBUG
    const arma::cx_mat44 S_S_ondiag_af_latticeBA = arma::kron(OneSiteSpinMatrices::get_S_z_af_latticeB(), OneSiteSpinMatrices::get_S_z_af_latticeA());
    assert(arma::norm(S_S_ondiag_af_latticeAB - S_S_ondiag_af_latticeBA) < 10 - 4);
#endif
    return S_S_ondiag_af_latticeAB;
}

arma::cx_mat44 TwoSitesSpinMatrices::get_S_S_offdiag_fm() {
    const arma::cx_mat44 S_S_ondiag_fm = 0.5 * (arma::kron(OneSiteSpinMatrices::get_S_plus_fm(), OneSiteSpinMatrices::get_S_minus_fm()) +
                                                arma::kron(OneSiteSpinMatrices::get_S_minus_fm(), OneSiteSpinMatrices::get_S_plus_fm()));
    return S_S_ondiag_fm;
}

arma::cx_mat44 TwoSitesSpinMatrices::get_S_S_offdiag_af() {
    const arma::cx_mat44 S_S_ondiagaf_latticeAB = 0.5 * (arma::kron(OneSiteSpinMatrices::get_S_plus_af_latticeA(), OneSiteSpinMatrices::get_S_minus_af_latticeB()) +
                                                         arma::kron(OneSiteSpinMatrices::get_S_minus_af_latticeA(), OneSiteSpinMatrices::get_S_plus_af_latticeB()));
#ifndef NDEBUG
    const arma::cx_mat44 S_S_ondiag_af_latticeBA = 0.5 * (arma::kron(OneSiteSpinMatrices::get_S_plus_af_latticeB(), OneSiteSpinMatrices::get_S_minus_af_latticeA()) +
                                                          arma::kron(OneSiteSpinMatrices::get_S_minus_af_latticeB(), OneSiteSpinMatrices::get_S_plus_af_latticeA()));
#endif
    return S_S_ondiagaf_latticeAB;
}

arma::cx_mat44 TwoSitesSpinMatrices::get_S_S_fm() {
    return get_S_S_ondiag_fm() + get_S_S_offdiag_fm();
}

arma::cx_mat44 TwoSitesSpinMatrices::get_S_S_af() {
    return get_S_S_ondiag_af() + get_S_S_offdiag_af();
}

}  // end of namespace monostar_hamiltonians

// #######################################################################
// ## SpinTwoSiteNamedMatrices                                          ##
// #######################################################################

namespace monostar_hamiltonians {

utility::Named<arma::cx_mat44> TwoSitesSpinNamedMatrices::get_S_S_ondiag_fm() {
    return utility::Named<arma::cx_mat44>(TwoSitesSpinMatrices::get_S_S_ondiag_fm(), "SᶻSᶻ");
}

utility::Named<arma::cx_mat44> TwoSitesSpinNamedMatrices::get_S_S_ondiag_af() {
    return utility::Named<arma::cx_mat44>(TwoSitesSpinMatrices::get_S_S_ondiag_af(), "SᶻSᶻ");
}

utility::Named<arma::cx_mat44> TwoSitesSpinNamedMatrices::get_S_S_offdiag_fm() {
    return utility::Named<arma::cx_mat44>(TwoSitesSpinMatrices::get_S_S_offdiag_fm(), "½(S⁺S⁻+S⁻S⁺)");
}

utility::Named<arma::cx_mat44> TwoSitesSpinNamedMatrices::get_S_S_offdiag_af() {
    return utility::Named<arma::cx_mat44>(TwoSitesSpinMatrices::get_S_S_offdiag_af(), "½(S⁺S⁻+S⁻S⁺)");
}

utility::Named<arma::cx_mat44> TwoSitesSpinNamedMatrices::get_S_S_fm() {
    return utility::Named<arma::cx_mat44>(TwoSitesSpinMatrices::get_S_S_fm(), "S⋅S");
}

utility::Named<arma::cx_mat44> TwoSitesSpinNamedMatrices::get_S_S_af() {
    return utility::Named<arma::cx_mat44>(TwoSitesSpinMatrices::get_S_S_af(), "S⋅S");
}

std::vector<utility::Named<arma::cx_mat44>> TwoSitesSpinNamedMatrices::matrices_for_average_calculations_fm() {
    return {TwoSitesSpinNamedMatrices::get_S_S_ondiag_fm(),
            TwoSitesSpinNamedMatrices::get_S_S_offdiag_fm(),
            TwoSitesSpinNamedMatrices::get_S_S_fm()};
}

std::vector<utility::Named<arma::cx_mat44>> TwoSitesSpinNamedMatrices::matrices_for_average_calculations_af() {
    return {TwoSitesSpinNamedMatrices::get_S_S_ondiag_af(),
            TwoSitesSpinNamedMatrices::get_S_S_offdiag_af(),
            TwoSitesSpinNamedMatrices::get_S_S_af()};
}

}  // end of namespace monostar_hamiltonians
