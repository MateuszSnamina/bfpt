#pragma once

#include<utility/named.hpp>

#include<armadillo>

#include<vector>

// #######################################################################
// ## OneSiteSpinMatrices                                               ##
// #######################################################################

namespace monostar_hamiltonians {

struct OneSiteSpinMatrices {
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

} // end of namespace monostar_hamiltonians

// #######################################################################
// ## OneSiteSpinNamedMatrices                                          ##
// #######################################################################

namespace monostar_hamiltonians {

struct OneSiteSpinNamedMatrices {
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

} // end of namespace monostar_hamiltonians

// #######################################################################
// ## TwoSitesSpinMatrices                                              ##
// #######################################################################

namespace monostar_hamiltonians {

struct TwoSitesSpinMatrices {
    static arma::cx_mat44 get_S_S_ondiag_fm();
    static arma::cx_mat44 get_S_S_ondiag_af();
    static arma::cx_mat44 get_S_S_offdiag_fm();
    static arma::cx_mat44 get_S_S_offdiag_af();
    static arma::cx_mat44 get_S_S_fm();
    static arma::cx_mat44 get_S_S_af();
};

} // end of namespace monostar_hamiltonians

// #######################################################################
// ## TwoSitesSpinNamedMatrices                                         ##
// #######################################################################

namespace monostar_hamiltonians {

struct TwoSitesSpinNamedMatrices {
    static utility::Named<arma::cx_mat44> get_S_S_ondiag_fm();
    static utility::Named<arma::cx_mat44> get_S_S_ondiag_af();
    static utility::Named<arma::cx_mat44> get_S_S_offdiag_fm();
    static utility::Named<arma::cx_mat44> get_S_S_offdiag_af();
    static utility::Named<arma::cx_mat44> get_S_S_fm();
    static utility::Named<arma::cx_mat44> get_S_S_af();
    static std::vector<utility::Named<arma::cx_mat44>> matrices_for_average_calculations_fm();
    static std::vector<utility::Named<arma::cx_mat44>> matrices_for_average_calculations_af();
};

} // end of namespace monostar_hamiltonians
