#include<model_monostar/monostar_hamiltonian.hpp>

// #######################################################################
// ## Helper function for preparing Hamiltonian12                       ##
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
        {{gs, gs}, {0.5 * J, {es, es}}}
    };
    return half_off_diag_info;
}

std::multimap<bfpt_common::StateKernel12<MonostarSiteState>, bfpt_common::CoupleInfoKernel12<MonostarSiteState>>
prepare_half_off_diag_info_for_fm(double J) {
    using RsultT = std::multimap<bfpt_common::StateKernel12<MonostarSiteState>, bfpt_common::CoupleInfoKernel12<MonostarSiteState>>;
    RsultT half_off_diag_info{
        {{gs, es}, {-0.5 * J, {es, gs}}}
    };
    return half_off_diag_info;
}

}

// #######################################################################
// ## prepare_hamiltonian_12_{af,fm}                                    ##
// #######################################################################

namespace model_monostar {

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_12_af(double J_classical, double J_quantum) {
    const auto diag_info = prepare_diag_info(J_classical);
    const auto half_off_diag_info = prepare_half_off_diag_info_for_af(J_quantum);
    return bfpt_common::HamiltonianKernel12<MonostarSiteState>{diag_info, half_off_diag_info};
}

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_12_fm(double J_classical, double J_quantum) {
    const auto diag_info = prepare_diag_info(J_classical);
    const auto half_off_diag_info = prepare_half_off_diag_info_for_fm(J_quantum);
    return bfpt_common::HamiltonianKernel12<MonostarSiteState>{diag_info, half_off_diag_info};
}

} // end of namespace model_monostar
