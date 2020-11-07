#pragma once

#include <utility/named.hpp>

#include <armadillo>

#include <vector>

// #######################################################################
// ## OneSiteOrbitalMatrices                                            ##
// #######################################################################

namespace so_hamiltonians {

struct OneSiteOrbitalMatrices {
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
// ## OneSiteOrbitalNamedMatrices                                       ##
// #######################################################################

namespace so_hamiltonians {

struct OneSiteOrbitalNamedMatrices {
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
    static std::vector<utility::Named<arma::cx_mat44>> matrices_for_average_calculations(double orbital_theta);
};

}  // end of namespace so_hamiltonians

// #######################################################################
// ## TwoSitesOrbitalMatrices                                           ##
// #######################################################################
/*
namespace so_hamiltonians {

struct TwoSitesOrbitalMatrices {
    static arma::cx_mat44 get_P_zz_in_zx_basis();
    static arma::cx_mat44 get_P_zx_sum_P_xz_in_zx_basis();
    static arma::cx_mat44 get_P_xx_in_zx_basis();
    static arma::cx_mat44 get_P_zz_in_ge_basis(double orbital_theta);
    static arma::cx_mat44 get_P_zx_sum_P_xz_in_ge_basis(double orbital_theta);
    static arma::cx_mat44 get_P_xx_in_ge_basis(double orbital_theta);
};

}  // end of namespace so_hamiltonians

// #######################################################################
// ## TwoSitesOrbitalNamedMatrices                                      ##
// #######################################################################

namespace so_hamiltonians {

struct TwoSitesOrbitalNamedMatrices {
    static utility::Named<arma::cx_mat44> get_P_zz_in_zx_basis();
    static utility::Named<arma::cx_mat44> get_P_zx_sum_P_xz_in_zx_basis();
    static utility::Named<arma::cx_mat44> get_P_xx_in_zx_basis();
    static utility::Named<arma::cx_mat44> get_P_zz_in_ge_basis(double orbital_theta);
    static utility::Named<arma::cx_mat44> get_P_zx_sum_P_xz_in_ge_basis(double orbital_theta);
    static utility::Named<arma::cx_mat44> get_P_xx_in_ge_basis(double orbital_theta);
    static std::vector<utility::Named<arma::cx_mat44>> matrices_for_average_calculations(double orbital_theta);
};

}  // end of namespace so_hamiltonians
*/
