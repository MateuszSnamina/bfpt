#ifndef BFPT_COMMON_DENSITY_OPERATOR_HPP
#define BFPT_COMMON_DENSITY_OPERATOR_HPP

#include <bfpt_common/operator_kernel.hpp>

#include <armadillo>

#include <map>
#include <complex>

// #######################################################################
// ## DensityOperator12                                                 ##
// #######################################################################

namespace bfpt_common {

template<typename SiteStateT>
using DensityOperator12 = std::map<std::pair<StateKernel12<SiteStateT>, StateKernel12<SiteStateT>>, std::complex<double>>;

template<typename SiteStateT>
arma::cx_mat
density_operator_12_to_matrix(
        const DensityOperator12<SiteStateT>& density_operator,
        std::vector<SiteStateT> ordered_site_states) {
    const arma::uword n_possible_site_states = ordered_site_states.size();
    arma::cx_mat rho(n_possible_site_states * n_possible_site_states, n_possible_site_states * n_possible_site_states);
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

} // end of namespace bfpt_common;

// #######################################################################
// ## DensityOperator1                                                  ##
// #######################################################################

namespace bfpt_common {

template<typename SiteState>
using DensityOperator1 = std::map<std::pair<StateKernel1<SiteState>, StateKernel1<SiteState>>, std::complex<double>>;


} // end of namespace bfpt_common;

#endif // CALCULATE_REDUCED_DENSITY_OPERATOR_HPP
