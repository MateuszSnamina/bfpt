#include<model_monostar/monostar_hamiltonian_kernel.hpp>

#include<cmath>

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
prepare_hamiltonian_kernel_12_fm(double J_classical, double J_quantum) {
    const auto diag_info = prepare_diag_info(J_classical);
    const auto half_off_diag_info = prepare_half_off_diag_info_for_fm(J_quantum);
    return bfpt_common::HamiltonianKernel12<MonostarSiteState>{diag_info, half_off_diag_info};
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

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_kernel_12_af_fm(double Pxx_coef, double Pxz_coef, double Pzz_coef, double theta_opt) {
    assert(Pxz_coef == 0.0); // TODO: implement nonzero case.
    using OnDiagInfoType = std::map<bfpt_common::StateKernel12<MonostarSiteState>, double>;
    using OffDiagInfoType = std::multimap<bfpt_common::StateKernel12<MonostarSiteState>, bfpt_common::CoupleInfoKernel12<MonostarSiteState>>;
    const double s2 = std::sin(theta_opt/2), c2 = std::cos(theta_opt/2);
    OnDiagInfoType on_diag_info {
        {{gs, gs}, +Pxx_coef * s2 * s2 * s2 * s2 + Pzz_coef * c2 * c2 * c2 * c2},
        {{gs, es}, +Pxx_coef * s2 * s2 * c2 * c2 + Pzz_coef * s2 * s2 * c2 * c2},
        {{es, gs}, +Pxx_coef * s2 * s2 * c2 * c2 + Pzz_coef * s2 * s2 * c2 * c2},
        {{es, es}, +Pxx_coef * c2 * c2 * c2 * c2 + Pzz_coef * s2 * s2 * s2 * s2}
    };
    OffDiagInfoType half_off_diag_info {
        {{gs, gs}, {+Pxx_coef * s2 * s2 * s2 * c2 - Pzz_coef * s2 * c2 * c2 * c2, {gs, es}}},
        {{gs, gs}, {+Pxx_coef * s2 * s2 * s2 * c2 - Pzz_coef * s2 * c2 * c2 * c2, {es, gs}}},
        {{gs, gs}, {+Pxx_coef * s2 * s2 * c2 * c2 + Pzz_coef * s2 * s2 * c2 * c2, {es, es}}},
        {{gs, es}, {+Pxx_coef * s2 * s2 * c2 * c2 + Pzz_coef * s2 * s2 * c2 * c2, {es, gs}}},
        {{gs, es}, {+Pxx_coef * s2 * c2 * c2 * c2 - Pzz_coef * s2 * s2 * s2 * c2, {es, es}}},
        {{es, gs}, {+Pxx_coef * s2 * c2 * c2 * c2 - Pzz_coef * s2 * s2 * s2 * c2, {es, es}}},
    };
    return bfpt_common::HamiltonianKernel12<MonostarSiteState>{on_diag_info, half_off_diag_info};
}

bfpt_common::HamiltonianKernel1<MonostarSiteState>
prepare_hamiltonian_kernel_1_fo(double Pdelta_coef, double theta_opt) {
    using OnDiagInfoType = std::map<bfpt_common::StateKernel1<MonostarSiteState>, double>;
    using OffDiagInfoType = std::multimap<bfpt_common::StateKernel1<MonostarSiteState>, bfpt_common::CoupleInfoKernel1<MonostarSiteState>>;
    const double s1 = std::sin(theta_opt), c1 = std::cos(theta_opt);
    OnDiagInfoType on_diag_info {
        {{gs}, +Pdelta_coef * s1},
        {{es}, -Pdelta_coef * s1},
    };
    OffDiagInfoType half_off_diag_info {
        {{gs}, {+Pdelta_coef * c1, {es}}},
    };
    return bfpt_common::HamiltonianKernel1<MonostarSiteState>{on_diag_info, half_off_diag_info};
}




} // end of namespace model_monostar
