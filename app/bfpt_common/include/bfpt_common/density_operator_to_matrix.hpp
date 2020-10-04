#ifndef BFPT_COMMON_DENSITY_OPERATOR_TO_MATRIX_HPP
#define BFPT_COMMON_DENSITY_OPERATOR_TO_MATRIX_HPP

#include <bfpt_common/density_operator.hpp>

#include <armadillo>

// #######################################################################
// ## density_operator_{1,12}_to_matrix                                 ##
// #######################################################################

namespace bfpt_common {

template<typename SiteStateT>
arma::cx_mat
density_operator_1_to_matrix(
        const DensityOperator12<SiteStateT>& density_operator,
        std::vector<SiteStateT> ordered_site_states) {
    const arma::uword n_possible_site_states = ordered_site_states.size();
    arma::cx_mat rho(n_possible_site_states, n_possible_site_states, arma::fill::zeros);
    for (arma::uword idx_bra_site_0 = 0; idx_bra_site_0 < n_possible_site_states; idx_bra_site_0++) {
        for (arma::uword idx_ket_site_0 = 0; idx_ket_site_0 < n_possible_site_states; idx_ket_site_0++) {
            const auto bra_site_0 = ordered_site_states[idx_bra_site_0];
            const auto ket_site_0 = ordered_site_states[idx_ket_site_0];
            const StateKernel12 bra_kenrel{bra_site_0};
            const StateKernel12 ket_kenrel{ket_site_0};
            const std::pair<StateKernel1<SiteStateT>, StateKernel1<SiteStateT>> density_matrix_indices{bra_kenrel, ket_kenrel};
            const std::complex<double> value = (density_operator.count(density_matrix_indices) ? density_operator.at(density_matrix_indices) : 0.0);
            rho(idx_bra_site_0, idx_ket_site_0) = value;
        }
    }
    return rho;
}

template<typename SiteStateT>
arma::cx_mat
density_operator_12_to_matrix(
        const DensityOperator12<SiteStateT>& density_operator,
        std::vector<SiteStateT> ordered_site_states) {
    const arma::uword n_possible_site_states = ordered_site_states.size();
    arma::cx_mat rho(n_possible_site_states * n_possible_site_states, n_possible_site_states * n_possible_site_states, arma::fill::zeros);
    for (arma::uword idx_bra_site_0 = 0; idx_bra_site_0 < n_possible_site_states; idx_bra_site_0++) {
        for (arma::uword idx_bra_site_1 = 0; idx_bra_site_1 < n_possible_site_states; idx_bra_site_1++) {
            for (arma::uword idx_ket_site_0 = 0; idx_ket_site_0 < n_possible_site_states; idx_ket_site_0++) {
                for (arma::uword idx_ket_site_1 = 0; idx_ket_site_1 < n_possible_site_states; idx_ket_site_1++) {
                    const auto bra_site_0 = ordered_site_states[idx_bra_site_0];
                    const auto bra_site_1 = ordered_site_states[idx_bra_site_1];
                    const auto ket_site_0 = ordered_site_states[idx_ket_site_0];
                    const auto ket_site_1 = ordered_site_states[idx_ket_site_1];
                    const StateKernel12 bra_kenrel{bra_site_0, bra_site_1};
                    const StateKernel12 ket_kenrel{ket_site_0, ket_site_1};
                    const std::pair<StateKernel12<SiteStateT>, StateKernel12<SiteStateT>> density_matrix_indices{bra_kenrel, ket_kenrel};
                    const std::complex<double> value = (density_operator.count(density_matrix_indices) ? density_operator.at(density_matrix_indices) : 0.0);
                    arma::uword idx_bra = n_possible_site_states * idx_bra_site_0 + idx_bra_site_1;
                    arma::uword idx_ket = n_possible_site_states * idx_ket_site_0 + idx_ket_site_1;
                    rho(idx_bra, idx_ket) = value;
                }
            }
        }
    }
    return rho;
}

} // end of namespace bfpt_common

#endif
