#pragma once

#include <kstate_trait/trait_site_state.hpp>

#include <type_traits>
#include <map>
#include <tuple>
#include <cassert>
#include <type_traits>

// #######################################################################
// ## OperatorKernel1                                                   ##
// #######################################################################

namespace chainkernel {

template <typename _SiteStateTraitT>
struct StateKernel1 {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    SiteStateT state_1;
};

template <typename SiteStateTraitT>
bool operator<(const StateKernel1<SiteStateTraitT>& lhs, const StateKernel1<SiteStateTraitT>& rhs) {
    return lhs.state_1 < rhs.state_1;
}

template <typename _SiteStateTraitT>
struct CoupleInfoKernel1 {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    StateKernel1<SiteStateTraitT> kernel_state;
    double coef;
};

template <typename _SiteStateTraitT>
class OperatorKernel1 {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);

   public:
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using OffDiagInfoT = std::multimap<StateKernel1<SiteStateTraitT>, CoupleInfoKernel1<SiteStateTraitT>>;
    using DiagInfoT = std::map<StateKernel1<SiteStateTraitT>, double>;

   public:
    OperatorKernel1() = default;
    OperatorKernel1(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info);
    const DiagInfoT _diag_info;
    const OffDiagInfoT _half_off_diag_info;
    const OffDiagInfoT _full_off_diag_info;
};

template <typename SiteStateTraitT>
std::multimap<StateKernel1<SiteStateTraitT>, CoupleInfoKernel1<SiteStateTraitT>>
half_off_diag_info_to_full_off_diag_info_1(
    const std::multimap<StateKernel1<SiteStateTraitT>, CoupleInfoKernel1<SiteStateTraitT>>& half_off_diag) {
    std::multimap<StateKernel1<SiteStateTraitT>, CoupleInfoKernel1<SiteStateTraitT>> full_off_diag_info;
    for (const auto node : half_off_diag) {
        const StateKernel1<SiteStateTraitT>& ket_1 = node.first;
        const CoupleInfoKernel1<SiteStateTraitT>& couple_info = node.second;
        const StateKernel1<SiteStateTraitT>& bra_1 = couple_info.kernel_state;
        const auto& kernel_coupling_coef = couple_info.coef;
        const auto& complementaty_ket_1 = bra_1;
        const auto& complementaty_bra_1 = ket_1;
        full_off_diag_info.insert({ket_1, {bra_1, kernel_coupling_coef}});
        full_off_diag_info.insert({complementaty_ket_1, {complementaty_bra_1, kernel_coupling_coef}});
    }
    return full_off_diag_info;
}

template <typename _SiteStateTraitT>
OperatorKernel1<_SiteStateTraitT>::OperatorKernel1(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info)
    : _diag_info(diag_info),
      _half_off_diag_info(half_off_diag_info),
      _full_off_diag_info(half_off_diag_info_to_full_off_diag_info_1(half_off_diag_info)) {
}

}  // namespace chainkernel

// #######################################################################
// ## OperatorKernel12                                                  ##
// #######################################################################

namespace chainkernel {

template <typename _SiteStateTraitT>
struct StateKernel12 {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    SiteStateT state_1;
    SiteStateT state_2;
};

template <typename SiteStateTraitT>
bool operator<(const StateKernel12<SiteStateTraitT>& lhs, const StateKernel12<SiteStateTraitT>& rhs) {
    return std::make_tuple(lhs.state_1, lhs.state_2) < std::make_tuple(rhs.state_1, rhs.state_2);
}

template <typename _SiteStateTraitT>
struct CoupleInfoKernel12 {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    StateKernel12<_SiteStateTraitT> kernel_state;
    double coef;
};

template <typename _SiteStateTraitT>
class OperatorKernel12 {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);

   public:
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using OffDiagInfoT = std::multimap<StateKernel12<SiteStateTraitT>, CoupleInfoKernel12<SiteStateTraitT>>;
    using DiagInfoT = std::map<StateKernel12<SiteStateTraitT>, double>;

   public:
    OperatorKernel12() = default;
    OperatorKernel12(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info);
    const DiagInfoT _diag_info;
    const OffDiagInfoT _half_off_diag_info;
    const OffDiagInfoT _full_off_diag_info;
};

template <typename SiteStateTraitT>
std::multimap<StateKernel12<SiteStateTraitT>, CoupleInfoKernel12<SiteStateTraitT>>
half_off_diag_info_to_full_off_diag_info_12(
    const std::multimap<StateKernel12<SiteStateTraitT>, CoupleInfoKernel12<SiteStateTraitT>>& half_off_diag) {
    std::multimap<StateKernel12<SiteStateTraitT>, CoupleInfoKernel12<SiteStateTraitT>> full_off_diag_info;
    for (const auto node : half_off_diag) {
        const StateKernel12<SiteStateTraitT>& ket_12 = node.first;
        const CoupleInfoKernel12<SiteStateTraitT>& couple_info = node.second;
        const StateKernel12<SiteStateTraitT>& bra_12 = couple_info.kernel_state;
        const auto& kernel_coupling_coef = couple_info.coef;
        const auto& complementaty_ket_12 = bra_12;
        const auto& complementaty_bra_12 = ket_12;
        full_off_diag_info.insert({ket_12, {bra_12, kernel_coupling_coef}});
        full_off_diag_info.insert({complementaty_ket_12, {complementaty_bra_12, kernel_coupling_coef}});
    }
    return full_off_diag_info;
}

template <typename _SiteStateTraitT>
OperatorKernel12<_SiteStateTraitT>::OperatorKernel12(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info)
    : _diag_info(diag_info),
      _half_off_diag_info(half_off_diag_info),
      _full_off_diag_info(half_off_diag_info_to_full_off_diag_info_12(half_off_diag_info)) {
}

}  // namespace chainkernel

// #######################################################################
// ## OperatorKernel123                                                 ##
// #######################################################################

namespace chainkernel {

template <typename _SiteStateTraitT>
struct StateKernel123 {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    SiteStateT state_1;
    SiteStateT state_2;
    SiteStateT state_3;
};

template <typename SiteStateTraitT>
bool operator<(const StateKernel123<SiteStateTraitT>& lhs, const StateKernel123<SiteStateTraitT>& rhs) {
    return std::make_tuple(lhs.state_1, lhs.state_2, lhs.state_3) < std::make_tuple(rhs.state_1, rhs.state_2, rhs.state_3);
}

template <typename _SiteStateTraitT>
struct CoupleInfoKernel123 {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    StateKernel123<_SiteStateTraitT> kernel_state;
    double coef;
};

template <typename _SiteStateTraitT>
class OperatorKernel123 {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);

   public:
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using OffDiagInfoT = std::multimap<StateKernel123<SiteStateTraitT>, CoupleInfoKernel123<SiteStateTraitT>>;
    using DiagInfoT = std::map<StateKernel123<SiteStateTraitT>, double>;

   public:
    OperatorKernel123() = default;
    OperatorKernel123(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info);
    const DiagInfoT _diag_info;
    const OffDiagInfoT _half_off_diag_info;
    const OffDiagInfoT _full_off_diag_info;
};

template <typename SiteStateTraitT>
std::multimap<StateKernel123<SiteStateTraitT>, CoupleInfoKernel123<SiteStateTraitT>>
half_off_diag_info_to_full_off_diag_info_123(
    const std::multimap<StateKernel123<SiteStateTraitT>, CoupleInfoKernel123<SiteStateTraitT>>& half_off_diag) {
    std::multimap<StateKernel123<SiteStateTraitT>, CoupleInfoKernel123<SiteStateTraitT>> full_off_diag_info;
    for (const auto node : half_off_diag) {
        const StateKernel123<SiteStateTraitT>& ket_123 = node.first;
        const CoupleInfoKernel123<SiteStateTraitT>& couple_info = node.second;
        const StateKernel123<SiteStateTraitT>& bra_123 = couple_info.kernel_state;
        const auto& kernel_coupling_coef = couple_info.coef;
        const auto& complementaty_ket_123 = bra_123;
        const auto& complementaty_bra_123 = ket_123;
        full_off_diag_info.insert({ket_123, {bra_123, kernel_coupling_coef}});
        full_off_diag_info.insert({complementaty_ket_123, {complementaty_bra_123, kernel_coupling_coef}});
    }
    return full_off_diag_info;
}

template <typename _SiteStateTraitT>
OperatorKernel123<_SiteStateTraitT>::OperatorKernel123(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info)
    : _diag_info(diag_info),
      _half_off_diag_info(half_off_diag_info),
      _full_off_diag_info(half_off_diag_info_to_full_off_diag_info_123(half_off_diag_info)) {
}

}  // namespace chainkernel

// #######################################################################
// ## OperatorKernel1234                                                ##
// #######################################################################

namespace chainkernel {

template <typename _SiteStateTraitT>
struct StateKernel1234 {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    SiteStateT state_1;
    SiteStateT state_2;
    SiteStateT state_3;
    SiteStateT state_4;
};

template <typename SiteStateTraitT>
bool operator<(const StateKernel1234<SiteStateTraitT>& lhs, const StateKernel1234<SiteStateTraitT>& rhs) {
    return std::make_tuple(lhs.state_1, lhs.state_2, lhs.state_3, lhs.state_4) < std::make_tuple(rhs.state_1, rhs.state_2, rhs.state_3, rhs.state_4);
}

template <typename _SiteStateTraitT>
struct CoupleInfoKernel1234 {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    StateKernel1234<_SiteStateTraitT> kernel_state;
    double coef;
};

template <typename _SiteStateTraitT>
class OperatorKernel1234 {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);

   public:
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using OffDiagInfoT = std::multimap<StateKernel1234<SiteStateTraitT>, CoupleInfoKernel1234<SiteStateTraitT>>;
    using DiagInfoT = std::map<StateKernel1234<SiteStateTraitT>, double>;

   public:
    OperatorKernel1234() = default;
    OperatorKernel1234(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info);
    const DiagInfoT _diag_info;
    const OffDiagInfoT _half_off_diag_info;
    const OffDiagInfoT _full_off_diag_info;
};

template <typename SiteStateTraitT>
std::multimap<StateKernel1234<SiteStateTraitT>, CoupleInfoKernel1234<SiteStateTraitT>>
half_off_diag_info_to_full_off_diag_info_1234(
    const std::multimap<StateKernel1234<SiteStateTraitT>, CoupleInfoKernel1234<SiteStateTraitT>>& half_off_diag) {
    std::multimap<StateKernel1234<SiteStateTraitT>, CoupleInfoKernel1234<SiteStateTraitT>> full_off_diag_info;
    for (const auto node : half_off_diag) {
        const StateKernel1234<SiteStateTraitT>& ket_1234 = node.first;
        const CoupleInfoKernel1234<SiteStateTraitT>& couple_info = node.second;
        const StateKernel1234<SiteStateTraitT>& bra_1234 = couple_info.kernel_state;
        const auto& kernel_coupling_coef = couple_info.coef;
        const auto& complementaty_ket_1234 = bra_1234;
        const auto& complementaty_bra_1234 = ket_1234;
        full_off_diag_info.insert({ket_1234, {bra_1234, kernel_coupling_coef}});
        full_off_diag_info.insert({complementaty_ket_1234, {complementaty_bra_1234, kernel_coupling_coef}});
    }
    return full_off_diag_info;
}

template <typename _SiteStateTraitT>
OperatorKernel1234<_SiteStateTraitT>::OperatorKernel1234(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info)
    : _diag_info(diag_info),
      _half_off_diag_info(half_off_diag_info),
      _full_off_diag_info(half_off_diag_info_to_full_off_diag_info_1234(half_off_diag_info)) {
}

}  // namespace chainkernel
