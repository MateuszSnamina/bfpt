#include<model_monostar/hamiltonian_kernel_fo.hpp>

#include<cmath>

// #######################################################################
// ## prepare_hamiltonian_kernel_{1,12}_fo                              ##
// #######################################################################

namespace model_monostar {

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_kernel_12_fo(double Pzz_coef, double Pxz_coef, double Pxx_coef, double orbital_theta) {
    using OnDiagInfoType = std::map<bfpt_common::StateKernel12<MonostarSiteState>, double>;
    using OffDiagInfoType = std::multimap<bfpt_common::StateKernel12<MonostarSiteState>, bfpt_common::CoupleInfoKernel12<MonostarSiteState>>;
    const double s2 = std::sin(orbital_theta/2);
    const double c2 = std::cos(orbital_theta/2);
    const double Pz_gg = +c2 * c2;
    const double Pz_ge = -s2 * c2;
    const double Pz_eg = -s2 * c2;
    const double Pz_ee = +s2 * s2;
    const double Px_gg = +s2 * s2;
    const double Px_ge = +s2 * c2;
    const double Px_eg = +s2 * c2;
    const double Px_ee = +c2 * c2;
    OnDiagInfoType on_diag_info {
        {{gs, gs}, + Pzz_coef * Pz_gg * Pz_gg + Pxz_coef * Px_gg * Pz_gg + Pxz_coef * Pz_gg * Px_gg + Pxx_coef * Px_gg * Px_gg},
        {{gs, es}, + Pzz_coef * Pz_gg * Pz_ee + Pxz_coef * Px_gg * Pz_ee + Pxz_coef * Pz_gg * Px_ee + Pxx_coef * Px_gg * Px_ee},
        {{es, gs}, + Pzz_coef * Pz_ee * Pz_gg + Pxz_coef * Px_ee * Pz_gg + Pxz_coef * Pz_ee * Px_gg + Pxx_coef * Px_ee * Px_gg},
        {{es, es}, + Pzz_coef * Pz_ee * Pz_ee + Pxz_coef * Px_ee * Pz_ee + Pxz_coef * Pz_ee * Px_ee + Pxx_coef * Px_ee * Px_ee}
    };
    OffDiagInfoType half_off_diag_info {
        {{gs, gs}, {{gs, es}, +Pzz_coef * Pz_gg * Pz_ge + Pxz_coef * Px_gg * Pz_ge + Pxz_coef * Pz_gg * Px_ge + Pxx_coef * Px_gg * Px_ge}},
        {{gs, gs}, {{es, gs}, +Pzz_coef * Pz_ge * Pz_gg + Pxz_coef * Px_ge * Pz_gg + Pxz_coef * Pz_ge * Px_gg + Pxx_coef * Px_ge * Px_gg}},
        {{gs, gs}, {{es, es}, +Pzz_coef * Pz_ge * Pz_ge + Pxz_coef * Px_ge * Pz_ge + Pxz_coef * Pz_ge * Px_ge + Pxx_coef * Px_ge * Px_ge}},
        {{gs, es}, {{es, gs}, +Pzz_coef * Pz_ge * Pz_eg + Pxz_coef * Px_ge * Pz_eg + Pxz_coef * Pz_ge * Px_eg + Pxx_coef * Px_ge * Px_eg}},
        {{gs, es}, {{es, es}, +Pzz_coef * Pz_ge * Pz_ee + Pxz_coef * Px_ge * Pz_ee + Pxz_coef * Pz_ge * Px_ee + Pxx_coef * Px_ge * Px_ee}},
        {{es, gs}, {{es, es}, +Pzz_coef * Pz_ee * Pz_ge + Pxz_coef * Px_ee * Pz_ge + Pxz_coef * Pz_ee * Px_ge + Pxx_coef * Px_ee * Px_ge}},
    };
    return bfpt_common::HamiltonianKernel12<MonostarSiteState>{on_diag_info, half_off_diag_info};
}

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_kernel_12_fo(const HamiltonianParamsFo& params, double orbital_theta) {
    return prepare_hamiltonian_kernel_12_fo(
                params.get_Pzz_coef(),
                params.get_Pxz_coef(),
                params.get_Pxx_coef(),
                orbital_theta);
}

bfpt_common::HamiltonianKernel1<MonostarSiteState>
prepare_hamiltonian_kernel_1_fo(double tau_z_coef, double tau_minus_coef, double orbital_theta) {
    using OnDiagInfoType = std::map<bfpt_common::StateKernel1<MonostarSiteState>, double>;
    using OffDiagInfoType = std::multimap<bfpt_common::StateKernel1<MonostarSiteState>, bfpt_common::CoupleInfoKernel1<MonostarSiteState>>;
    const double s1 = std::sin(orbital_theta), c1 = std::cos(orbital_theta);
    const double tau_z_gg = +c1;
    const double tau_z_ge = -s1;
    const double tau_z_ee = -c1;
    const double tau_minus_gg = -s1;
    const double tau_minus_ge = -c1;
    const double tau_minus_ee = +s1;
    OnDiagInfoType on_diag_info {
        {{gs}, tau_z_coef * tau_z_gg + tau_minus_coef * tau_minus_gg},
        {{es}, tau_z_coef * tau_z_ee + tau_minus_coef * tau_minus_ee},
    };
    OffDiagInfoType half_off_diag_info {
        {{gs}, {{es}, tau_z_coef * tau_z_ge + tau_minus_coef * tau_minus_ge}},
    };
    return bfpt_common::HamiltonianKernel1<MonostarSiteState>{on_diag_info, half_off_diag_info};
}

bfpt_common::HamiltonianKernel1<MonostarSiteState>
prepare_hamiltonian_kernel_1_fo(const HamiltonianParamsFo& params, double orbital_theta) {
    return prepare_hamiltonian_kernel_1_fo(
                params.get_tau_z_coef(),
                params.get_tau_minus_coef(),
                orbital_theta);
}

} // end of namespace model_monostar
