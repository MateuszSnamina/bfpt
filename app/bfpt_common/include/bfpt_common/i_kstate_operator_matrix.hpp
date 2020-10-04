#ifndef BFPT_COMMON_I_KSTATE_OPERATOR_HPP
#define BFPT_COMMON_I_KSTATE_OPERATOR_HPP

#include <kstate/basis.hpp>

#include <armadillo>

// #######################################################################
// ## IKstateOperatorMatrix                                             ##
// #######################################################################

namespace bfpt_common {

/*
 * This inteface is to providing kstate operator matrix builder.
 * Technically, the operator is build from matrix `G`
 * such as `H = G + G^\dagger`.
 *
 * The class provides a member function
 * for filling a single column in `G` matrix.
 * (The arg name `kn_operator_builder_matrix` is for the matrix G).
 *
 */

template<typename _KstateTraitT>
class IKstateOperatorMatrix {
    static_assert(kstate::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);
public:
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename _KstateTraitT::KstateT;
    using SiteStateTraitT = typename KstateT::SiteStateTraitT;
    using SiteStateT = typename KstateT::SiteStateT;
    using BasisT = kstate::Basis<KstateTraitT>;
public:
    virtual void fill_kn_operator_builder_matrix_coll(
            const BasisT& basis,
            size_t n_col,
            arma::sp_cx_mat& kn_operator_builder_matrix,
            const unsigned k_n) const = 0;
};

}  // namespace bfpt_common

#endif
