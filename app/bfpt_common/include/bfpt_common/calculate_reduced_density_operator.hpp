#ifndef BFPT_COMMON_CALCULATE_REDUCED_DENSITY_OPERATOR_HPP
#define BFPT_COMMON_CALCULATE_REDUCED_DENSITY_OPERATOR_HPP

#include <bfpt_common/density_operator.hpp>

#include <kstate/basis.hpp>
#include <kstate/trait_kstate.hpp>

#include <boost/range/adaptor/sliced.hpp>
#include <boost/range/algorithm/equal.hpp>

#include <armadillo>
#include <map>
#include <complex>

// #######################################################################
// ## DensityOperator12                                                 ##
// #######################################################################

namespace bfpt_common {

template<typename KstateTraitT>
DensityOperator12<typename KstateTraitT::KstateT::SiteStateTraitT>
calculate_reduced_density_operator_12(
        kstate::Basis<KstateTraitT>& basis,
        const arma::cx_vec& eigen_vector) {
    static_assert(kstate::IsTraitKstate<KstateTraitT>::value);
    static_assert(KstateTraitT::is_kstate_trait);
    using KstateT = typename KstateTraitT::KstateT;
    using SiteStateTraitT = typename KstateT::SiteStateTraitT;
    //using SiteStateT = typename KstateT::SiteStateT;
    DensityOperator12<SiteStateTraitT> result;
    const auto n_sites = basis.n_sites();
    assert(eigen_vector.n_rows == basis.size());
    for (arma::uword bra_kstate_idx = 0; bra_kstate_idx < eigen_vector.n_rows; bra_kstate_idx++) {
        for (arma::uword ket_kstate_idx = 0; ket_kstate_idx < eigen_vector.n_rows; ket_kstate_idx++) {
            const auto bra_kstate_ptr = basis.vec_index()[bra_kstate_idx];
            const auto ket_kstate_ptr = basis.vec_index()[ket_kstate_idx];
            assert(bra_kstate_ptr);
            assert(ket_kstate_ptr);
            const auto bra_kstate_range = KstateTraitT::to_range(*bra_kstate_ptr);
            const auto ket_kstate_range = KstateTraitT::to_range(*ket_kstate_ptr);
            assert(boost::size(bra_kstate_range) == n_sites);
            assert(boost::size(ket_kstate_range) == n_sites);
            for (unsigned bra_shift = 0; bra_shift < n_sites; bra_shift++) {
                for (unsigned ket_shift = 0; ket_shift < n_sites; ket_shift++) {
                    const auto bra_kstate_range_rotated = bra_kstate_range | extension::boost::adaptors::rotated(bra_shift);
                    const auto ket_kstate_range_rotated = ket_kstate_range | extension::boost::adaptors::rotated(ket_shift);
                    const auto bra_kstate_range_rotated_rest = bra_kstate_range_rotated | boost::adaptors::sliced(2, n_sites);
                    const auto ket_kstate_range_rotated_rest = ket_kstate_range_rotated | boost::adaptors::sliced(2, n_sites);
                    if (boost::equal(bra_kstate_range_rotated_rest, ket_kstate_range_rotated_rest)) {
                        const auto bra_site_0 = *std::next(std::begin(bra_kstate_range_rotated), 0);
                        const auto bra_site_1 = *std::next(std::begin(bra_kstate_range_rotated), 1);
                        const auto ket_site_0 = *std::next(std::begin(ket_kstate_range_rotated), 0);
                        const auto ket_site_1 = *std::next(std::begin(ket_kstate_range_rotated), 1);
                        const StateKernel12<SiteStateTraitT> bra_kenrel{bra_site_0, bra_site_1};
                        const StateKernel12<SiteStateTraitT> ket_kenrel{ket_site_0, ket_site_1};
                        const std::pair<StateKernel12<SiteStateTraitT>, StateKernel12<SiteStateTraitT>> density_matrix_indices{bra_kenrel, ket_kenrel};
                        const std::complex<double> value = std::conj(eigen_vector(bra_kstate_idx)) * eigen_vector(ket_kstate_idx);
                        const double pre_norm = bra_kstate_ptr->norm_factor() * ket_kstate_ptr->norm_factor();
                        if (result.count(density_matrix_indices)) {
                            result[density_matrix_indices] += pre_norm * value;
                        } else {
                            result[density_matrix_indices] = pre_norm * value;
                        }
                    }
                } // end of ket_shift loop
            } // end of bra_shift loop
        } // end of ket_kstate_idx loop
    } // end of bra_kstate_idx loop
    return result;
}

} // end of namespace bfpt_common;

// #######################################################################
// ## DensityOperator1                                                  ##
// #######################################################################

namespace bfpt_common {

template<typename KstateTraitT>
DensityOperator1<typename KstateTraitT::KstateT::SiteStateTraitT>
calculate_reduced_density_operator_1(
        kstate::Basis<KstateTraitT>& basis,
        const arma::cx_vec& eigen_vector) {
    static_assert(kstate::IsTraitKstate<KstateTraitT>::value);
    static_assert(KstateTraitT::is_kstate_trait);
    using KstateT = typename KstateTraitT::KstateT;
    using SiteStateTraitT = typename KstateT::SiteStateTraitT;
    //using SiteStateT = typename KstateT::SiteStateT;
    DensityOperator1<SiteStateTraitT> result;
    const auto n_sites = basis.n_sites();
    assert(eigen_vector.n_rows == basis.size());
    for (arma::uword bra_kstate_idx = 0; bra_kstate_idx < eigen_vector.n_rows; bra_kstate_idx++) {
        for (arma::uword ket_kstate_idx = 0; ket_kstate_idx < eigen_vector.n_rows; ket_kstate_idx++) {
            const auto bra_kstate_ptr = basis.vec_index()[bra_kstate_idx];
            const auto ket_kstate_ptr = basis.vec_index()[ket_kstate_idx];
            assert(bra_kstate_ptr);
            assert(ket_kstate_ptr);
            const auto bra_kstate_range = (*bra_kstate_ptr).to_range();
            const auto ket_kstate_range = (*ket_kstate_ptr).to_range();
            assert(boost::size(bra_kstate_range) == n_sites);
            assert(boost::size(ket_kstate_range) == n_sites);
            for (unsigned bra_shift = 0; bra_shift < n_sites; bra_shift++) {
                for (unsigned ket_shift = 0; ket_shift < n_sites; ket_shift++) {
                    const auto bra_kstate_range_rotated = bra_kstate_range | extension::boost::adaptors::rotated(bra_shift);
                    const auto ket_kstate_range_rotated = ket_kstate_range | extension::boost::adaptors::rotated(ket_shift);
                    const auto bra_kstate_range_rotated_rest = bra_kstate_range_rotated | boost::adaptors::sliced(1, n_sites);
                    const auto ket_kstate_range_rotated_rest = ket_kstate_range_rotated | boost::adaptors::sliced(1, n_sites);
                    if (boost::equal(bra_kstate_range_rotated_rest, ket_kstate_range_rotated_rest)) {
                        const auto bra_site_0 = *std::next(std::begin(bra_kstate_range_rotated), 0);
                        const auto ket_site_0 = *std::next(std::begin(ket_kstate_range_rotated), 0);
                        const StateKernel1<SiteStateTraitT> bra_kenrel{bra_site_0};
                        const StateKernel1<SiteStateTraitT> ket_kenrel{ket_site_0};
                        const std::pair<StateKernel1<SiteStateTraitT>, StateKernel1<SiteStateTraitT>> density_matrix_indices{bra_kenrel, ket_kenrel};
                        const std::complex<double> value = std::conj(eigen_vector(bra_kstate_idx)) * eigen_vector(ket_kstate_idx);
                        const double pre_norm = bra_kstate_ptr->norm_factor() * ket_kstate_ptr->norm_factor();
                        if (result.count(density_matrix_indices)) {
                            result[density_matrix_indices] += pre_norm * value;
                        } else {
                            result[density_matrix_indices] = pre_norm * value;
                        }
                    }
                } // end of ket_shift loop
            } // end of bra_shift loop
        } // end of ket_kstate_idx loop
    } // end of bra_kstate_idx loop
    return result;
}

} // end of namespace bfpt_common;

#endif // CALCULATE_REDUCED_DENSITY_OPERATOR_HPP
