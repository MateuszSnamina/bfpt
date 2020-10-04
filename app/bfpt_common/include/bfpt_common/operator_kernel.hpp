#ifndef BFPT_COMMON_HAMILTONIAN_KERNEL_HPP
#define BFPT_COMMON_HAMILTONIAN_KERNEL_HPP

#include<type_traits>
#include<map>
#include<tuple>
#include<cassert>
#include<type_traits>

// #######################################################################
// ## OperatorKernel1                                                   ##
// #######################################################################

namespace bfpt_common {

template<typename _SiteStateTrait>
struct StateKernel1 {
    static_assert(_SiteStateTrait::is_site_state_trait);
    using SiteStateTrait = _SiteStateTrait;
    using SiteStateT = typename _SiteStateTrait::SiteStateT;
    SiteStateT state_1;
};

template<typename SiteStateTrait>
bool operator<(const StateKernel1<SiteStateTrait>& lhs, const StateKernel1<SiteStateTrait>& rhs) {
    return lhs.state_1 < rhs.state_1;
}

template<typename _SiteStateTrait>
struct CoupleInfoKernel1 {
    static_assert(_SiteStateTrait::is_site_state_trait);
    using SiteStateTrait = _SiteStateTrait;
    using SiteStateT = typename _SiteStateTrait::SiteStateT;
    StateKernel1<SiteStateTrait> kernel_state;
    double coef;
};

template<typename _SiteStateTrait>
class OperatorKernel1 {
    static_assert(_SiteStateTrait::is_site_state_trait);
public:
    using SiteStateTrait = _SiteStateTrait;
    using SiteStateT = typename _SiteStateTrait::SiteStateT;
    using OffDiagInfoT = std::multimap<StateKernel1<SiteStateTrait>, CoupleInfoKernel1<SiteStateTrait>>;
    using DiagInfoT = std::map<StateKernel1<SiteStateTrait>, double> ;
public:
    OperatorKernel1(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info);
    const DiagInfoT _diag_info;
    const OffDiagInfoT _half_off_diag_info;
    const OffDiagInfoT _full_off_diag_info;
};

template<typename SiteStateTrait>
std::multimap<StateKernel1<SiteStateTrait>, CoupleInfoKernel1<SiteStateTrait>>
half_off_diag_info_to_full_off_diag_info_1(
        const std::multimap<StateKernel1<SiteStateTrait>, CoupleInfoKernel1<SiteStateTrait>>& half_off_diag) {
    std::multimap<StateKernel1<SiteStateTrait>, CoupleInfoKernel1<SiteStateTrait>> full_off_diag_info;
    for (const auto node : half_off_diag) {
        const StateKernel1<SiteStateTrait>& ket_1 = node.first;
        const CoupleInfoKernel1<SiteStateTrait>& couple_info = node.second;
        const StateKernel1<SiteStateTrait>& bra_1 = couple_info.kernel_state;
        const auto& kernel_coupling_coef = couple_info.coef;
        const auto& complementaty_ket_1 = bra_1;
        const auto& complementaty_bra_1 = ket_1;
        full_off_diag_info.insert({ket_1, {bra_1, kernel_coupling_coef, }});
        full_off_diag_info.insert({complementaty_ket_1, {complementaty_bra_1, kernel_coupling_coef}});
    }
    return full_off_diag_info;
}

template<typename _SiteStateTrait>
OperatorKernel1<_SiteStateTrait>::OperatorKernel1(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info) :
    _diag_info(diag_info),
    _half_off_diag_info(half_off_diag_info),
    _full_off_diag_info(half_off_diag_info_to_full_off_diag_info_1(half_off_diag_info)) {
}

}

// #######################################################################
// ## OperatorKernel12                                                  ##
// #######################################################################

namespace bfpt_common {

template<typename _SiteStateTrait>
struct StateKernel12 {
    static_assert(_SiteStateTrait::is_site_state_trait);
    using SiteStateTrait = _SiteStateTrait;
    using SiteStateT = typename _SiteStateTrait::SiteStateT;
    SiteStateT state_1;
    SiteStateT state_2;
};

template<typename SiteStateTrait>
bool operator<(const StateKernel12<SiteStateTrait>& lhs, const StateKernel12<SiteStateTrait>& rhs) {
    return std::make_tuple(lhs.state_1, lhs.state_2) < std::make_tuple(rhs.state_1, rhs.state_2);
}

template<typename _SiteStateTrait>
struct CoupleInfoKernel12 {
    static_assert(_SiteStateTrait::is_site_state_trait);
    using SiteStateTrait = _SiteStateTrait;
    using SiteStateT = typename _SiteStateTrait::SiteStateT;
    StateKernel12<_SiteStateTrait> kernel_state;
    double coef;
};

template<typename _SiteStateTrait>
class OperatorKernel12 {
    static_assert(_SiteStateTrait::is_site_state_trait);
public:
    using SiteStateTrait = _SiteStateTrait;
    using SiteStateT = typename _SiteStateTrait::SiteStateT;
    using OffDiagInfoT = std::multimap<StateKernel12<SiteStateTrait>, CoupleInfoKernel12<SiteStateTrait>>;
    using DiagInfoT = std::map<StateKernel12<SiteStateTrait>, double> ;
public:
    OperatorKernel12(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info);
    const DiagInfoT _diag_info;
    const OffDiagInfoT _half_off_diag_info;
    const OffDiagInfoT _full_off_diag_info;
};

template<typename SiteStateTrait>
std::multimap<StateKernel12<SiteStateTrait>, CoupleInfoKernel12<SiteStateTrait>>
half_off_diag_info_to_full_off_diag_info_12(
        const std::multimap<StateKernel12<SiteStateTrait>, CoupleInfoKernel12<SiteStateTrait>>& half_off_diag) {
    std::multimap<StateKernel12<SiteStateTrait>, CoupleInfoKernel12<SiteStateTrait>> full_off_diag_info;
    for (const auto node : half_off_diag) {
        const StateKernel12<SiteStateTrait>& ket_12 = node.first;
        const CoupleInfoKernel12<SiteStateTrait>& couple_info = node.second;
        const StateKernel12<SiteStateTrait>& bra_12 = couple_info.kernel_state;
        const auto& kernel_coupling_coef = couple_info.coef;
        const auto& complementaty_ket_12 = bra_12;
        const auto& complementaty_bra_12 = ket_12;
        full_off_diag_info.insert({ket_12, {bra_12, kernel_coupling_coef}});
        full_off_diag_info.insert({complementaty_ket_12, {complementaty_bra_12, kernel_coupling_coef}});
    }
    return full_off_diag_info;
}

template<typename _SiteStateTrait>
OperatorKernel12<_SiteStateTrait>::OperatorKernel12(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info) :
    _diag_info(diag_info),
    _half_off_diag_info(half_off_diag_info),
    _full_off_diag_info(half_off_diag_info_to_full_off_diag_info_12(half_off_diag_info)) {
}

}

#endif
