#pragma once

#include<armadillo>

struct OrbitalSiteMatrices {

    static arma::cx_mat22 get_P_z_in_zx_basis();
    static arma::cx_mat22 get_P_x_in_zx_basis();
    static arma::cx_mat22 get_P_plus_in_zx_basis();
    static arma::cx_mat22 get_P_minus_in_zx_basis();
    static arma::cx_mat22 get_tau_z_in_zx_basis();
    static arma::cx_mat22 get_tau_x_in_zx_basis();
    static arma::cx_mat22 get_tau_plus_in_zx_basis();
    static arma::cx_mat22 get_tau_minus_in_zx_basis();

    static arma::cx_mat22 get_beta_from_zx_to_ge(double orbital_theta);

    static arma::cx_mat22 get_P_z_in_ge_basis(double orbital_theta);
    static arma::cx_mat22 get_P_x_in_ge_basis(double orbital_theta);
    static arma::cx_mat22 get_P_plus_in_ge_basis(double orbital_theta);
    static arma::cx_mat22 get_P_minus_in_ge_basis(double orbital_theta);
    static arma::cx_mat22 get_tau_z_in_ge_basis(double orbital_theta);
    static arma::cx_mat22 get_tau_x_in_ge_basis(double orbital_theta);
    static arma::cx_mat22 get_tau_plus_in_ge_basis(double orbital_theta);
    static arma::cx_mat22 get_tau_minus_in_ge_basis(double orbital_theta);

};
