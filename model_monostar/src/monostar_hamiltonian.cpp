#include<model_monostar/monostar_hamiltonian.hpp>

// #######################################################################
// ## Helper function for preparing Hamiltonian12                       ##
// #######################################################################

namespace model_monostar {

std::multimap<bfpt_common::SiteStatePair<MonostarSiteState>, bfpt_common::CoupleInfo<MonostarSiteState>>
prepare_half_off_diag_info_for_af(double J) {
    using RsultT = std::multimap<bfpt_common::SiteStatePair<MonostarSiteState>, bfpt_common::CoupleInfo<MonostarSiteState>>;
    RsultT half_off_diag_info{
        {{gs, gs}, {0.5 * J, {es,es}}}
    };
    return half_off_diag_info;
}

std::map<bfpt_common::SiteStatePair<MonostarSiteState>, double>
prepare_diag_info(double J) {
    using RsultT = std::map<bfpt_common::SiteStatePair<MonostarSiteState>, double>;
    RsultT diag_info{
        {{gs, gs}, -J * 0.25},
        {{gs, es}, +J * 0.25},
        {{es, gs}, +J * 0.25},
        {{es, es}, -J * 0.25}
    };
    return diag_info;
}

bfpt_common::Hamiltonian12<MonostarSiteState>
prepare_hamiltonian_12(double J_classical, double J_quantum) {
    const auto diag_info = prepare_diag_info(J_classical);
    const auto half_off_diag_info = prepare_half_off_diag_info_for_af(J_quantum);
    return bfpt_common::Hamiltonian12<MonostarSiteState>{diag_info, half_off_diag_info};
}

} // end of namespace model_monostar
