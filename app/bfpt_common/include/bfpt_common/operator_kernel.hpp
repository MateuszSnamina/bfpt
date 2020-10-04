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

template<typename _SiteStateT>
struct StateKernel1 {
    static_assert(!std::is_array_v<_SiteStateT>);
    static_assert(!std::is_function_v<_SiteStateT>);
    static_assert(!std::is_void_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_null_pointer_v<std::decay<_SiteStateT>>);
    static_assert(std::is_enum_v<std::decay<_SiteStateT>> || std::is_union_v<std::decay<_SiteStateT>> || std::is_class_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_pointer_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_member_object_pointer_v<_SiteStateT>);
    static_assert(!std::is_member_function_pointer_v<_SiteStateT>);
    static_assert(!std::is_const_v<_SiteStateT>);
    static_assert(!std::is_volatile_v<_SiteStateT>);
    static_assert(!std::is_reference_v<_SiteStateT>);
    using SiteStateT = _SiteStateT;
    SiteStateT state_1;
};

template<typename SiteStateT>
bool operator<(const StateKernel1<SiteStateT>& lhs, const StateKernel1<SiteStateT>& rhs) {
    return lhs.state_1 < rhs.state_1;
}

template<typename _SiteStateT>
struct CoupleInfoKernel1 {
    static_assert(!std::is_array_v<_SiteStateT>);
    static_assert(!std::is_function_v<_SiteStateT>);
    static_assert(!std::is_void_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_null_pointer_v<std::decay<_SiteStateT>>);
    static_assert(std::is_enum_v<std::decay<_SiteStateT>> || std::is_union_v<std::decay<_SiteStateT>> || std::is_class_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_pointer_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_member_object_pointer_v<_SiteStateT>);
    static_assert(!std::is_member_function_pointer_v<_SiteStateT>);
    static_assert(!std::is_const_v<_SiteStateT>);
    static_assert(!std::is_volatile_v<_SiteStateT>);
    static_assert(!std::is_reference_v<_SiteStateT>);
    using SiteStateT = _SiteStateT;
    StateKernel1<SiteStateT> kernel_state;
    double coef;
};

template<typename _SiteStateT>
class OperatorKernel1 {
    static_assert(!std::is_array_v<_SiteStateT>);
    static_assert(!std::is_function_v<_SiteStateT>);
    static_assert(!std::is_void_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_null_pointer_v<std::decay<_SiteStateT>>);
    static_assert(std::is_enum_v<std::decay<_SiteStateT>> || std::is_union_v<std::decay<_SiteStateT>> || std::is_class_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_pointer_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_member_object_pointer_v<_SiteStateT>);
    static_assert(!std::is_member_function_pointer_v<_SiteStateT>);
    static_assert(!std::is_const_v<_SiteStateT>);
    static_assert(!std::is_volatile_v<_SiteStateT>);
    static_assert(!std::is_reference_v<_SiteStateT>);
public:
    using SiteStateT = _SiteStateT;
    using OffDiagInfoT = std::multimap<StateKernel1<SiteStateT>, CoupleInfoKernel1<SiteStateT>>;
    using DiagInfoT = std::map<StateKernel1<SiteStateT>, double> ;
public:
    OperatorKernel1(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info);
    const DiagInfoT _diag_info;
    const OffDiagInfoT _half_off_diag_info;
    const OffDiagInfoT _full_off_diag_info;
};

template<typename SiteState>
std::multimap<StateKernel1<SiteState>, CoupleInfoKernel1<SiteState>>
half_off_diag_info_to_full_off_diag_info_1(
        const std::multimap<StateKernel1<SiteState>, CoupleInfoKernel1<SiteState>>& half_off_diag) {
    std::multimap<StateKernel1<SiteState>, CoupleInfoKernel1<SiteState>> full_off_diag_info;
    for (const auto node : half_off_diag) {
        const StateKernel1<SiteState>& ket_1 = node.first;
        const CoupleInfoKernel1<SiteState>& couple_info = node.second;
        const StateKernel1<SiteState>& bra_1 = couple_info.kernel_state;
        const auto& kernel_coupling_coef = couple_info.coef;
        const auto& complementaty_ket_1 = bra_1;
        const auto& complementaty_bra_1 = ket_1;
        full_off_diag_info.insert({ket_1, {bra_1, kernel_coupling_coef, }});
        full_off_diag_info.insert({complementaty_ket_1, {complementaty_bra_1, kernel_coupling_coef}});
    }
    return full_off_diag_info;
}

template<typename _SiteState>
OperatorKernel1<_SiteState>::OperatorKernel1(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info) :
    _diag_info(diag_info),
    _half_off_diag_info(half_off_diag_info),
    _full_off_diag_info(half_off_diag_info_to_full_off_diag_info_1(half_off_diag_info)) {
}

}

// #######################################################################
// ## OperatorKernel12                                                  ##
// #######################################################################

namespace bfpt_common {

template<typename _SiteStateT>
struct StateKernel12 {
    static_assert(!std::is_array_v<_SiteStateT>);
    static_assert(!std::is_function_v<_SiteStateT>);
    static_assert(!std::is_void_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_null_pointer_v<std::decay<_SiteStateT>>);
    static_assert(std::is_enum_v<std::decay<_SiteStateT>> || std::is_union_v<std::decay<_SiteStateT>> || std::is_class_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_pointer_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_member_object_pointer_v<_SiteStateT>);
    static_assert(!std::is_member_function_pointer_v<_SiteStateT>);
    static_assert(!std::is_const_v<_SiteStateT>);
    static_assert(!std::is_volatile_v<_SiteStateT>);
    static_assert(!std::is_reference_v<_SiteStateT>);
    using SiteStateT = _SiteStateT;
    SiteStateT state_1;
    SiteStateT state_2;
};

template<typename SiteStateT>
bool operator<(const StateKernel12<SiteStateT>& lhs, const StateKernel12<SiteStateT>& rhs) {
    return std::make_tuple(lhs.state_1, lhs.state_2) < std::make_tuple(rhs.state_1, rhs.state_2);
}

template<typename _SiteStateT>
struct CoupleInfoKernel12 {
    static_assert(!std::is_array_v<_SiteStateT>);
    static_assert(!std::is_function_v<_SiteStateT>);
    static_assert(!std::is_void_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_null_pointer_v<std::decay<_SiteStateT>>);
    static_assert(std::is_enum_v<std::decay<_SiteStateT>> || std::is_union_v<std::decay<_SiteStateT>> || std::is_class_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_pointer_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_member_object_pointer_v<_SiteStateT>);
    static_assert(!std::is_member_function_pointer_v<_SiteStateT>);
    static_assert(!std::is_const_v<_SiteStateT>);
    static_assert(!std::is_volatile_v<_SiteStateT>);
    static_assert(!std::is_reference_v<_SiteStateT>);
    using SiteStateT = _SiteStateT;
    StateKernel12<SiteStateT> kernel_state;
    double coef;
};

template<typename _SiteStateT>
class OperatorKernel12 {
    static_assert(!std::is_array_v<_SiteStateT>);
    static_assert(!std::is_function_v<_SiteStateT>);
    static_assert(!std::is_void_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_null_pointer_v<std::decay<_SiteStateT>>);
    static_assert(std::is_enum_v<std::decay<_SiteStateT>> || std::is_union_v<std::decay<_SiteStateT>> || std::is_class_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_pointer_v<std::decay<_SiteStateT>>);
    static_assert(!std::is_member_object_pointer_v<_SiteStateT>);
    static_assert(!std::is_member_function_pointer_v<_SiteStateT>);
    static_assert(!std::is_const_v<_SiteStateT>);
    static_assert(!std::is_volatile_v<_SiteStateT>);
    static_assert(!std::is_reference_v<_SiteStateT>);
public:
    using SiteStateT = _SiteStateT;
    using OffDiagInfoT = std::multimap<StateKernel12<SiteStateT>, CoupleInfoKernel12<SiteStateT>>;
    using DiagInfoT = std::map<StateKernel12<SiteStateT>, double> ;
public:
    OperatorKernel12(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info);
    const DiagInfoT _diag_info;
    const OffDiagInfoT _half_off_diag_info;
    const OffDiagInfoT _full_off_diag_info;
};

template<typename SiteState>
std::multimap<StateKernel12<SiteState>, CoupleInfoKernel12<SiteState>>
half_off_diag_info_to_full_off_diag_info_12(
        const std::multimap<StateKernel12<SiteState>, CoupleInfoKernel12<SiteState>>& half_off_diag) {
    std::multimap<StateKernel12<SiteState>, CoupleInfoKernel12<SiteState>> full_off_diag_info;
    for (const auto node : half_off_diag) {
        const StateKernel12<SiteState>& ket_12 = node.first;
        const CoupleInfoKernel12<SiteState>& couple_info = node.second;
        const StateKernel12<SiteState>& bra_12 = couple_info.kernel_state;
        const auto& kernel_coupling_coef = couple_info.coef;
        const auto& complementaty_ket_12 = bra_12;
        const auto& complementaty_bra_12 = ket_12;
        full_off_diag_info.insert({ket_12, {bra_12, kernel_coupling_coef}});
        full_off_diag_info.insert({complementaty_ket_12, {complementaty_bra_12, kernel_coupling_coef}});
    }
    return full_off_diag_info;
}

template<typename _SiteState>
OperatorKernel12<_SiteState>::OperatorKernel12(DiagInfoT diag_info, OffDiagInfoT half_off_diag_info) :
    _diag_info(diag_info),
    _half_off_diag_info(half_off_diag_info),
    _full_off_diag_info(half_off_diag_info_to_full_off_diag_info_12(half_off_diag_info)) {
}

}

#endif
