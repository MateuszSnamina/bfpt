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
    static_assert(!std::is_const<_KstateT>::value);
    static_assert(!std::is_volatile<_KstateT>::value);
    static_assert(!std::is_reference<_KstateT>::value);
    static_assert(kstate::is_base_of_template_v<_KstateT, kstate::Kstate>);
public:
    using KstateT = _KstateT;
    using SiteStateT = typename kstate::remove_cvref_t<KstateT>::SiteType;
    using BasisT = kstate::Basis<KstateT>;

public:
    virtual arma::sp_cx_mat make_kn_hamiltonian_matrix(
            const BasisT& basis,
            const unsigned k_n) const = 0;
};

}  // namespace bfpt_common

#endif
