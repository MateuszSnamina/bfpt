#ifndef BFPT_COMMON_I_DYNAMIC_UNIQUE_KSTATE_HAMILTONIAN_HPP
#define BFPT_COMMON_I_DYNAMIC_UNIQUE_KSTATE_HAMILTONIAN_HPP

#include <kstate/basis.hpp>
#include <kstate/kstate_concrete.hpp>

#include <armadillo>

// #######################################################################
// ## IDynamicUniqueKstateHamiltonian                                   ##
// #######################################################################

namespace bfpt_common {

template <typename _SiteType>
class IDynamicUniqueKstateHamiltonian {
   public:
    virtual arma::sp_cx_mat make_kn_hamiltonian_matrix(
        const kstate::Basis<kstate::DynamicUniqueKstate<_SiteType>>& basis,
        const unsigned k_n) const = 0;
};

}  // namespace bfpt_common

#endif
