#ifndef BFPT_COMMON_I_KSTATE_OPERATOR_HPP
#define BFPT_COMMON_I_KSTATE_OPERATOR_HPP

#include <kstate/remove_cvref.hpp>
#include <kstate/kstate_abstract.hpp>
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

template<typename _KstateT>
class IKstateOperatorMatrix {
    static_assert(!std::is_array_v<_KstateT>);
    static_assert(!std::is_function_v<_KstateT>);
    static_assert(!std::is_void_v<std::decay<_KstateT>>);
    static_assert(!std::is_null_pointer_v<std::decay<_KstateT>>);
    static_assert(!std::is_enum_v<std::decay<_KstateT>>);
    static_assert(!std::is_union_v<std::decay<_KstateT>>);
    static_assert(std::is_class_v<std::decay<_KstateT>>);
    static_assert(!std::is_pointer_v<std::decay<_KstateT>>);
    static_assert(!std::is_member_object_pointer_v<_KstateT>);
    static_assert(!std::is_member_function_pointer_v<_KstateT>);
    static_assert(!std::is_const_v<_KstateT>);
    static_assert(!std::is_volatile_v<_KstateT>);
    static_assert(!std::is_reference_v<_KstateT>);
    static_assert(kstate::is_base_of_template_v<_KstateT, kstate::Kstate>);
public:
    using KstateT = _KstateT;
    using SiteStateT = typename kstate::remove_cvref_t<KstateT>::SiteType;
    using BasisT = kstate::Basis<KstateT>;

public:
    virtual void fill_kn_operator_builder_matrix_coll(
            const BasisT& basis,
            size_t n_col,
            arma::sp_cx_mat& kn_operator_builder_matrix,
            const unsigned k_n) const = 0;
};

}  // namespace bfpt_common

#endif
