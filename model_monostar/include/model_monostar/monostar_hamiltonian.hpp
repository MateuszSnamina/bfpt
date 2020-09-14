#ifndef MODEL_MONOSTAR_MONOSTAR_HAMILTONIAN_HPP
#define MODEL_MONOSTAR_MONOSTAR_HAMILTONIAN_HPP

#include <model_monostar/monostar_basis.hpp>
#include <model_monostar/monostar_kstate.hpp>

#include <bfpt_common/i_dynamic_unique_kstate_hamiltonian.hpp>
#include <bfpt_common/i_dynamic_unique_kstate_populator.hpp>

#include <armadillo>

#include <type_traits>
#include <map>
#include <tuple>

// #######################################################################
// ## DynamicMonostarHamiltonian                                        ##
// #######################################################################

namespace model_monostar {

template<typename _SiteStateT>
struct SiteStatePair {
    using SiteStateT = _SiteStateT;
    static_assert(!std::is_reference_v<SiteStateT>, "SiteState must not be a reference type.");
    static_assert(!std::is_const_v<SiteStateT>, "SiteState must not be a cosnt type.");
    static_assert(!std::is_volatile_v<SiteStateT>, "SiteState must not be a volatile type.");
    SiteStateT state_1;
    SiteStateT state_2;
};

template<typename SiteStateT>
bool operator<(const SiteStatePair<SiteStateT>& lhs, const SiteStatePair<SiteStateT>& rhs) {
    return std::make_tuple(lhs.state_1, lhs.state_2) < std::make_tuple(rhs.state_1, rhs.state_2);
}

template<typename _SiteStateT>
struct CoupleInfo {
    using SiteStateT = _SiteStateT;
    double coef;
    SiteStatePair<SiteStateT> state12;
};

template<typename _SiteStateT>
class Hamiltonian12 {
public:
    using SiteStateT = _SiteStateT;
    using OffDiagInfoT = std::multimap<SiteStatePair<SiteStateT>, CoupleInfo<SiteStateT>>;
    using DiagInfoT = std::map<SiteStatePair<SiteStateT>, double> ;
public:
    Hamiltonian12(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info);
    const DiagInfoT _diag_info;
    const OffDiagInfoT _half_off_diag_info;
    const OffDiagInfoT _full_off_diag_info;
};

template<typename SiteState>
std::multimap<SiteStatePair<SiteState>, CoupleInfo<SiteState>>
half_off_diag_info_to_full_off_diag_info(
        const std::multimap<SiteStatePair<SiteState>, CoupleInfo<SiteState>>& half_off_diag) {
    std::multimap<SiteStatePair<SiteState>, CoupleInfo<SiteState>> full_off_diag_info;
    for (const auto node : half_off_diag) {
        const SiteStatePair<SiteState>& ket_12 = node.first;
        const CoupleInfo<SiteState>& couple_info = node.second;
        const SiteStatePair<SiteState>& bra_12 = couple_info.state12;
        const auto& kernel_coupling_coef = couple_info.coef;
        const auto& complementaty_ket_12 = bra_12;
        const auto& complementaty_bra_12 = ket_12;
        full_off_diag_info.insert({ket_12, {kernel_coupling_coef, bra_12}});
        full_off_diag_info.insert({complementaty_ket_12, {kernel_coupling_coef, complementaty_bra_12}});
    }
    return full_off_diag_info;
}

template<typename _SiteState>
Hamiltonian12<_SiteState>::Hamiltonian12(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info) :
    _diag_info(diag_info),
    _half_off_diag_info(half_off_diag_info),
    _full_off_diag_info(half_off_diag_info_to_full_off_diag_info(half_off_diag_info)) {
}

// #######################################################################
// ## DynamicMonostarHamiltonian                                        ##
// #######################################################################

class DynamicMonostarHamiltonian : public bfpt_common::IDynamicUniqueKstatePopulator<MonostarSiteState>,
        public bfpt_common::IDynamicUniqueKstateHamiltonian<MonostarSiteState> {
public:
    using SiteState = MonostarSiteState;
public:
    DynamicMonostarHamiltonian(const size_t n_sites, Hamiltonian12<SiteState> hamiltonian_12);
    // Generates all conjugated states:
    void push_back_coupled_states_to_basis(
            const DynamicMonostarUniqueKstate& generator,
            DynamicMonostarUniqueKstateBasis& basis) const override;
    // Generates hamiltonian matrix:
    arma::sp_cx_mat make_kn_hamiltonian_matrix(
            const DynamicMonostarUniqueKstateBasis& basis,
            const unsigned k_n) const override;

private:
    void fill_kn_hamiltonian_matrix_coll(
            const DynamicMonostarUniqueKstateBasis& basis,
            size_t n_col,
            arma::sp_cx_mat& kn_hamiltonian_matrix,
            const unsigned k_n) const;
    void fill_kn_hamiltonian_matrix(
            const DynamicMonostarUniqueKstateBasis& basis,
            arma::sp_cx_mat& kn_hamiltonian_matrix,
            const unsigned k_n) const;

private:
    const size_t _n_sites;
    const Hamiltonian12<SiteState> _hamiltonian_12;
};

}  // namespace model_monostar

#endif
