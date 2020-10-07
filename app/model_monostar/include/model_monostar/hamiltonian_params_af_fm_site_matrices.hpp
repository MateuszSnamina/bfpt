#pragma once

#include<utility/named.hpp>

#include<armadillo>

#include<vector>

// #######################################################################
// ## SpinSiteMatrices                                                  ##
// #######################################################################

struct SpinSiteMatrices {
    static arma::cx_mat22 get_S_z_fm();
    static arma::cx_mat22 get_S_z_af_latticeA();
    static arma::cx_mat22 get_S_z_af_latticeB();
    static arma::cx_mat22 get_S_plus_fm();
    static arma::cx_mat22 get_S_plus_af_latticeA();
    static arma::cx_mat22 get_S_plus_af_latticeB();
    static arma::cx_mat22 get_S_minus_fm();
    static arma::cx_mat22 get_S_minus_af_latticeA();
    static arma::cx_mat22 get_S_minus_af_latticeB();
};

// #######################################################################
// ## SpinSiteNamedMatrices                                             ##
// #######################################################################

struct SpinSiteNamedMatrices {
    static utility::Named<arma::cx_mat22> get_S_z_fm();
    static utility::Named<arma::cx_mat22> get_S_z_af_latticeA();
    static utility::Named<arma::cx_mat22> get_S_z_af_latticeB();
    static utility::Named<arma::cx_mat22> get_S_plus_fm();
    static utility::Named<arma::cx_mat22> get_S_plus_af_latticeA();
    static utility::Named<arma::cx_mat22> get_S_plus_af_latticeB();
    static utility::Named<arma::cx_mat22> get_S_minus_fm();
    static utility::Named<arma::cx_mat22> get_S_minus_af_latticeA();
    static utility::Named<arma::cx_mat22> get_S_minus_af_latticeB();
    static std::vector<utility::Named<arma::cx_mat22>> site_matrices_for_average_calculations_af_fm();
};

// #######################################################################
// ## SpinTwoSiteMatrices                                               ##
// #######################################################################

struct SpinTwoSiteMatrices {
    static arma::cx_mat44 get_S_S_ondiag_fm();
    static arma::cx_mat44 get_S_S_ondiag_af();
    static arma::cx_mat44 get_S_S_offdiag_fm();
    static arma::cx_mat44 get_S_S_offdiag_af();
    static arma::cx_mat44 get_S_S_fm();
    static arma::cx_mat44 get_S_S_af();
};

// #######################################################################
// ## SpinTwoSiteNamedMatrices                                          ##
// #######################################################################

struct SpinTwoSiteNamedMatrices {
    static utility::Named<arma::cx_mat44> get_S_S_ondiag_fm();
    static utility::Named<arma::cx_mat44> get_S_S_ondiag_af();
    static utility::Named<arma::cx_mat44> get_S_S_offdiag_fm();
    static utility::Named<arma::cx_mat44> get_S_S_offdiag_af();
    static utility::Named<arma::cx_mat44> get_S_S_fm();
    static utility::Named<arma::cx_mat44> get_S_S_af();
    static std::vector<utility::Named<arma::cx_mat44>> two_site_matrices_for_average_calculations_fm();
    static std::vector<utility::Named<arma::cx_mat44>> two_site_matrices_for_average_calculations_af();
};
