#ifndef BFPT_COMMON_I_DYNAMIC_UNIQUE_KSTATE_HAMILTONIAN_HPP
#define BFPT_COMMON_I_DYNAMIC_UNIQUE_KSTATE_HAMILTONIAN_HPP

#include <kstate/remove_cvref.hpp>
#include <kstate/kstate_abstract.hpp>
#include <kstate/basis.hpp>

#include <armadillo>

// #######################################################################
// ## IDynamicUniqueKstateHamiltonian                                   ##
// #######################################################################

namespace bfpt_common {

template<typename _KstateT>
class IKstateHamiltonian {
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
    virtual arma::sp_cx_mat make_kn_hamiltonian_matrix(
            const BasisT& basis,
            const unsigned k_n,
            unsigned n_threads) const = 0;
};

}  // namespace bfpt_common

#endif
