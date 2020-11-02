#include <monostar_hamiltonians/hamiltonian_kernel_af_fm.hpp>

#include <cmath>

// #######################################################################
// ## helper functions...                                               ##
// #######################################################################

namespace {

using namespace monostar_system;

std::map<chainkernel::StateKernel12<MonostarSiteStateTrait>, double>
prepare_diag_info(double J) {
    using RsultT = std::map<chainkernel::StateKernel12<MonostarSiteStateTrait>, double>;
    RsultT diag_info{
        {{gs, gs}, -J * 0.25},
        {{gs, es}, +J * 0.25},
        {{es, gs}, +J * 0.25},
        {{es, es}, -J * 0.25}};
    return diag_info;
}

std::multimap<
    chainkernel::StateKernel12<MonostarSiteStateTrait>,
    chainkernel::CoupleInfoKernel12<MonostarSiteStateTrait>>
prepare_half_off_diag_info_for_af(double J) {
    using RsultT = std::multimap<
        chainkernel::StateKernel12<MonostarSiteStateTrait>,
        chainkernel::CoupleInfoKernel12<MonostarSiteStateTrait>>;
    RsultT half_off_diag_info{
        {{gs, gs}, {{es, es}, 0.5 * J}}};
    return half_off_diag_info;
}

std::multimap<chainkernel::StateKernel12<MonostarSiteStateTrait>, chainkernel::CoupleInfoKernel12<MonostarSiteStateTrait>>
prepare_half_off_diag_info_for_fm(double J) {
    using RsultT = std::multimap<chainkernel::StateKernel12<MonostarSiteStateTrait>, chainkernel::CoupleInfoKernel12<MonostarSiteStateTrait>>;
    RsultT half_off_diag_info{
        {{gs, es}, {{es, gs}, -0.5 * J}}};
    return half_off_diag_info;
}

}  // namespace

// #######################################################################
// ## prepare_hamiltonian_kernel_{1,12}_{af,fm}                         ##
// #######################################################################

namespace monostar_hamiltonians {

chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_af(double J_classical, double J_quantum) {
    const auto diag_info = prepare_diag_info(J_classical);
    const auto half_off_diag_info = prepare_half_off_diag_info_for_af(J_quantum);
    return chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>{diag_info, half_off_diag_info};
}

chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_af(const HamiltonianParamsAfFm& params) {
    return prepare_hamiltonian_kernel_12_af(params.get_J_classical(), params.get_J_quantum());
}

chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_fm(double J_classical, double J_quantum) {
    const auto diag_info = prepare_diag_info(J_classical);
    const auto half_off_diag_info = prepare_half_off_diag_info_for_fm(J_quantum);
    return chainkernel::OperatorKernel12<MonostarSiteStateTrait>{diag_info, half_off_diag_info};
}

chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_fm(const HamiltonianParamsAfFm& params) {
    return prepare_hamiltonian_kernel_12_fm(params.get_J_classical(), params.get_J_quantum());
}

chainkernel::OperatorKernel1<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1_af_fm(double B) {
    using namespace monostar_system;
    using OnDiagInfoType = std::map<chainkernel::StateKernel1<MonostarSiteStateTrait>, double>;
    using OffDiagInfoType = std::multimap<
        chainkernel::StateKernel1<MonostarSiteStateTrait>,
        chainkernel::CoupleInfoKernel1<MonostarSiteStateTrait>>;
    OnDiagInfoType on_diag_info{
        {{gs}, -B * 0.5},
        {{es}, +B * 0.5},
    };
    OffDiagInfoType half_off_diag_info{};
    return chainkernel::OperatorKernel1<MonostarSiteStateTrait>{on_diag_info, half_off_diag_info};
}

chainkernel::OperatorKernel1<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1_af_fm(const HamiltonianParamsAfFm& params) {
    return prepare_hamiltonian_kernel_1_af_fm(params.get_B());
}

}  // end of namespace monostar_hamiltonians
