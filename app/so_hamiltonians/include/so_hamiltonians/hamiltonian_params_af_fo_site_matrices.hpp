#pragma once

#include <utility/named.hpp>

#include <armadillo>

#include <vector>

// #######################################################################
// ## OneSiteSpinOrbitalMatrices                                        ##
// #######################################################################

namespace so_hamiltonians {

struct OneSiteSpinOrbitalMatrices {
    // spin operators:
    static arma::cx_mat44 get_S_z_af_latticeA();
    static arma::cx_mat44 get_S_z_af_latticeB();
    static arma::cx_mat44 get_S_plus_af_latticeA();
    static arma::cx_mat44 get_S_plus_af_latticeB();
    static arma::cx_mat44 get_S_minus_af_latticeA();
    static arma::cx_mat44 get_S_minus_af_latticeB();
    // orbital operators:
    static arma::cx_mat44 get_P_z_in_zx_basis();
    static arma::cx_mat44 get_P_x_in_zx_basis();
    static arma::cx_mat44 get_P_plus_in_zx_basis();
    static arma::cx_mat44 get_P_minus_in_zx_basis();
    static arma::cx_mat44 get_tau_z_in_zx_basis();
    static arma::cx_mat44 get_tau_x_in_zx_basis();
    static arma::cx_mat44 get_tau_plus_in_zx_basis();
    static arma::cx_mat44 get_tau_minus_in_zx_basis();
    static arma::cx_mat44 get_beta_from_zx_to_ge(double orbital_theta);
    static arma::cx_mat44 get_P_z_in_ge_basis(double orbital_theta);
    static arma::cx_mat44 get_P_x_in_ge_basis(double orbital_theta);
    static arma::cx_mat44 get_P_plus_in_ge_basis(double orbital_theta);
    static arma::cx_mat44 get_P_minus_in_ge_basis(double orbital_theta);
    static arma::cx_mat44 get_tau_z_in_ge_basis(double orbital_theta);
    static arma::cx_mat44 get_tau_x_in_ge_basis(double orbital_theta);
    static arma::cx_mat44 get_tau_plus_in_ge_basis(double orbital_theta);
    static arma::cx_mat44 get_tau_minus_in_ge_basis(double orbital_theta);
};

}  // end of namespace so_hamiltonians

// #######################################################################
// ## OneSiteSpinOrbitalNamedMatrices                                   ##
// #######################################################################

namespace so_hamiltonians {

struct OneSiteSpinOrbitalNamedMatrices {
    // spin operators:
    static utility::Named<arma::cx_mat44> get_S_z_af_latticeA();
    static utility::Named<arma::cx_mat44> get_S_z_af_latticeB();
    static utility::Named<arma::cx_mat44> get_S_plus_af_latticeA();
    static utility::Named<arma::cx_mat44> get_S_plus_af_latticeB();
    static utility::Named<arma::cx_mat44> get_S_minus_af_latticeA();
    static utility::Named<arma::cx_mat44> get_S_minus_af_latticeB();
    // orbital operators:
    static utility::Named<arma::cx_mat44> get_P_z_in_zx_basis();
    static utility::Named<arma::cx_mat44> get_P_x_in_zx_basis();
    static utility::Named<arma::cx_mat44> get_P_plus_in_zx_basis();
    static utility::Named<arma::cx_mat44> get_P_minus_in_zx_basis();
    static utility::Named<arma::cx_mat44> get_tau_z_in_zx_basis();
    static utility::Named<arma::cx_mat44> get_tau_x_in_zx_basis();
    static utility::Named<arma::cx_mat44> get_tau_plus_in_zx_basis();
    static utility::Named<arma::cx_mat44> get_tau_minus_in_zx_basis();
    static utility::Named<arma::cx_mat44> get_P_z_in_ge_basis(double orbital_theta);
    static utility::Named<arma::cx_mat44> get_P_x_in_ge_basis(double orbital_theta);
    static utility::Named<arma::cx_mat44> get_P_plus_in_ge_basis(double orbital_theta);
    static utility::Named<arma::cx_mat44> get_P_minus_in_ge_basis(double orbital_theta);
    static utility::Named<arma::cx_mat44> get_tau_z_in_ge_basis(double orbital_theta);
    static utility::Named<arma::cx_mat44> get_tau_x_in_ge_basis(double orbital_theta);
    static utility::Named<arma::cx_mat44> get_tau_plus_in_ge_basis(double orbital_theta);
    static utility::Named<arma::cx_mat44> get_tau_minus_in_ge_basis(double orbital_theta);
    // for average calculations:
    static std::vector<utility::Named<arma::cx_mat44>> matrices_for_average_calculations(double orbital_theta);
};

}  // end of namespace so_hamiltonians

// #######################################################################
// ## TwoSitesSpinOrbitalMatrices                                       ##
// #######################################################################


namespace arma::extension {
    using cx_mat1616 = ::arma::cx_mat::fixed<16, 16>;
}


namespace so_hamiltonians {

struct TwoSitesSpinOrbitalMatrices {
    // spin operators:
    static arma::extension::cx_mat1616 get_S_S_ondiag_fm();
    static arma::extension::cx_mat1616 get_S_S_ondiag_af();
    static arma::extension::cx_mat1616 get_S_S_offdiag_fm();
    static arma::extension::cx_mat1616 get_S_S_offdiag_af();
    static arma::extension::cx_mat1616 get_S_S_fm();
    static arma::extension::cx_mat1616 get_S_S_af();
    // orbital operators:
    static arma::extension::cx_mat1616 get_P_zz_in_zx_basis();
    static arma::extension::cx_mat1616 get_P_zx_sum_P_xz_in_zx_basis();
    static arma::extension::cx_mat1616 get_P_xx_in_zx_basis();
    static arma::extension::cx_mat1616 get_P_zz_in_ge_basis(double orbital_theta);
    static arma::extension::cx_mat1616 get_P_zx_sum_P_xz_in_ge_basis(double orbital_theta);
    static arma::extension::cx_mat1616 get_P_xx_in_ge_basis(double orbital_theta);
};

}  // end of namespace so_hamiltonians

// #######################################################################
// ## TwoSitesSpinOrbitalNamedMatrices                                  ##
// #######################################################################

namespace so_hamiltonians {

struct TwoSitesSpinOrbitalNamedMatrices {
    // spin operators:
    static utility::Named<arma::extension::cx_mat1616> get_S_S_ondiag_fm();
    static utility::Named<arma::extension::cx_mat1616> get_S_S_ondiag_af();
    static utility::Named<arma::extension::cx_mat1616> get_S_S_offdiag_fm();
    static utility::Named<arma::extension::cx_mat1616> get_S_S_offdiag_af();
    static utility::Named<arma::extension::cx_mat1616> get_S_S_fm();
    static utility::Named<arma::extension::cx_mat1616> get_S_S_af();
    // orbital operators:
    static utility::Named<arma::extension::cx_mat1616> get_P_zz_in_zx_basis();
    static utility::Named<arma::extension::cx_mat1616> get_P_zx_sum_P_xz_in_zx_basis();
    static utility::Named<arma::extension::cx_mat1616> get_P_xx_in_zx_basis();
    static utility::Named<arma::extension::cx_mat1616> get_P_zz_in_ge_basis(double orbital_theta);
    static utility::Named<arma::extension::cx_mat1616> get_P_zx_sum_P_xz_in_ge_basis(double orbital_theta);
    static utility::Named<arma::extension::cx_mat1616> get_P_xx_in_ge_basis(double orbital_theta);
    // for average calculations:
    static std::vector<utility::Named<arma::extension::cx_mat1616>> matrices_for_average_calculations(double orbital_theta);
};

}  // end of namespace so_hamiltonians
