#ifndef BFPT_COMMON_HAMILTONIAN_12_HPP
#define BFPT_COMMON_HAMILTONIAN_12_HPP

#include<type_traits>
#include<map>
#include<tuple>
#include<cassert>

// #######################################################################
// ## Hamiltonian12                                                     ##
// #######################################################################

namespace bfpt_common {

template<typename _SiteStateT>
struct SiteStatePair {
    static_assert(!std::is_reference_v<_SiteStateT>, "SiteState must not be a reference type.");
    static_assert(!std::is_const_v<_SiteStateT>, "SiteState must not be a cosnt type.");
    static_assert(!std::is_volatile_v<_SiteStateT>, "SiteState must not be a volatile type.");
    using SiteStateT = _SiteStateT;
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

}

#endif
