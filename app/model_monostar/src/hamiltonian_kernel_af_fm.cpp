#include<model_monostar/hamiltonian_kernel_af_fm.hpp>

#include<cmath>

// #######################################################################
// ## helper functions...                                               ##
// #######################################################################

namespace {

using namespace model_monostar;

std::map<bfpt_common::StateKernel12<MonostarSiteState>, double>
prepare_diag_info(double J) {
    using RsultT = std::map<bfpt_common::StateKernel12<MonostarSiteState>, double>;
    RsultT diag_info{
        {{gs, gs}, -J * 0.25},
        {{gs, es}, +J * 0.25},
        {{es, gs}, +J * 0.25},
        {{es, es}, -J * 0.25}
    };
    return diag_info;
}

std::multimap<bfpt_common::StateKernel12<MonostarSiteState>, bfpt_common::CoupleInfoKernel12<MonostarSiteState>>
prepare_half_off_diag_info_for_af(double J) {
    using RsultT = std::multimap<bfpt_common::StateKernel12<MonostarSiteState>, bfpt_common::CoupleInfoKernel12<MonostarSiteState>>;
    RsultT half_off_diag_info{
        {{gs, gs}, {{es, es}, 0.5 * J}}
    };
    return half_off_diag_info;
}

std::multimap<bfpt_common::StateKernel12<MonostarSiteState>, bfpt_common::CoupleInfoKernel12<MonostarSiteState>>
prepare_half_off_diag_info_for_fm(double J) {
    using RsultT = std::multimap<bfpt_common::StateKernel12<MonostarSiteState>, bfpt_common::CoupleInfoKernel12<MonostarSiteState>>;
    RsultT half_off_diag_info{
        {{gs, es}, {{es, gs}, -0.5 * J}}
    };
    return half_off_diag_info;
}

}

// #######################################################################
// ## prepare_hamiltonian_kernel_{1,12}_{af,fm}                         ##
// #######################################################################

namespace model_monostar {

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_kernel_12_af(double J_classical, double J_quantum) {
    const auto diag_info = prepare_diag_info(J_classical);
    const auto half_off_diag_info = prepare_half_off_diag_info_for_af(J_quantum);
    return bfpt_common::HamiltonianKernel12<MonostarSiteState>{diag_info, half_off_diag_info};
}

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_kernel_12_af(const HamiltonianAfFmParams& params) {
    return prepare_hamiltonian_kernel_12_af(params.get_J_classical(), params.get_J_quantum());
}

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_kernel_12_fm(double J_classical, double J_quantum) {
    const auto diag_info = prepare_diag_info(J_classical);
    const auto half_off_diag_info = prepare_half_off_diag_info_for_fm(J_quantum);
    return bfpt_common::HamiltonianKernel12<MonostarSiteState>{diag_info, half_off_diag_info};
}

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_kernel_12_fm(const HamiltonianAfFmParams& params) {
    return prepare_hamiltonian_kernel_12_fm(params.get_J_classical(), params.get_J_quantum());
}

bfpt_common::HamiltonianKernel1<MonostarSiteState>
prepare_hamiltonian_kernel_1_af_fm(double B) {
    using OnDiagInfoType = std::map<bfpt_common::StateKernel1<MonostarSiteState>, double>;
    using OffDiagInfoType = std::multimap<bfpt_common::StateKernel1<MonostarSiteState>, bfpt_common::CoupleInfoKernel1<MonostarSiteState>>;
    OnDiagInfoType on_diag_info{
        {{gs}, -B * 0.5},
        {{es}, +B * 0.5},
    };
    OffDiagInfoType half_off_diag_info{
    };
    return bfpt_common::HamiltonianKernel1<MonostarSiteState>{on_diag_info, half_off_diag_info};
}

bfpt_common::HamiltonianKernel1<MonostarSiteState>
prepare_hamiltonian_kernel_1_af_fm(const HamiltonianAfFmParams& params) {
    return prepare_hamiltonian_kernel_1_af_fm(params.get_B());
}

} // end of namespace model_monostar
