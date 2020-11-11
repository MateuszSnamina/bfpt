#include <so_hamiltonians/hamiltonian_kernel_af_fo.hpp>

#include <cmath>

// #######################################################################
// ## prepare_hamiltonian_kernel_{1,12}_fo                              ##
// #######################################################################

namespace so_hamiltonians {

chainkernel::OperatorKernel12<so_system::SoSiteStateTrait>
prepare_hamiltonian_kernel_12_fo(
    double ss_coef,
    double Pzz_coef, double Pxz_coef, double Pxx_coef,
    double ss_Pzz_coef, double ss_Pxz_coef, double ss_Pxx_coef,
    double orbital_theta) {
    using namespace so_system;
    using OnDiagInfoType = std::map<chainkernel::StateKernel12<SoSiteStateTrait>, double>;
    using OffDiagInfoType = std::multimap<chainkernel::StateKernel12<SoSiteStateTrait>, chainkernel::CoupleInfoKernel12<SoSiteStateTrait>>;
    const double s2 = std::sin(orbital_theta / 2);
    const double c2 = std::cos(orbital_theta / 2);
    const double Pz_gg = +c2 * c2;
    const double Pz_ge = -s2 * c2;
    const double Pz_eg = -s2 * c2;
    const double Pz_ee = +s2 * s2;
    const double Px_gg = +s2 * s2;
    const double Px_ge = +s2 * c2;
    const double Px_eg = +s2 * c2;
    const double Px_ee = +c2 * c2;
    // On-diag:
    const double orbit_gg_gg = +Pzz_coef * Pz_gg * Pz_gg + Pxz_coef * Px_gg * Pz_gg + Pxz_coef * Pz_gg * Px_gg + Pxx_coef * Px_gg * Px_gg;
    const double orbit_ge_ge = +Pzz_coef * Pz_gg * Pz_ee + Pxz_coef * Px_gg * Pz_ee + Pxz_coef * Pz_gg * Px_ee + Pxx_coef * Px_gg * Px_ee;
    const double orbit_eg_eg = +Pzz_coef * Pz_ee * Pz_gg + Pxz_coef * Px_ee * Pz_gg + Pxz_coef * Pz_ee * Px_gg + Pxx_coef * Px_ee * Px_gg;
    const double orbit_ee_ee = +Pzz_coef * Pz_ee * Pz_ee + Pxz_coef * Px_ee * Pz_ee + Pxz_coef * Pz_ee * Px_ee + Pxx_coef * Px_ee * Px_ee;
    const double spinorbit_gg_gg = +ss_Pzz_coef * Pz_gg * Pz_gg + ss_Pxz_coef * Px_gg * Pz_gg + ss_Pxz_coef * Pz_gg * Px_gg + ss_Pxx_coef * Px_gg * Px_gg;
    const double spinorbit_ge_ge = +ss_Pzz_coef * Pz_gg * Pz_ee + ss_Pxz_coef * Px_gg * Pz_ee + ss_Pxz_coef * Pz_gg * Px_ee + ss_Pxx_coef * Px_gg * Px_ee;
    const double spinorbit_eg_eg = +ss_Pzz_coef * Pz_ee * Pz_gg + ss_Pxz_coef * Px_ee * Pz_gg + ss_Pxz_coef * Pz_ee * Px_gg + ss_Pxx_coef * Px_ee * Px_gg;
    const double spinorbit_ee_ee = +ss_Pzz_coef * Pz_ee * Pz_ee + ss_Pxz_coef * Px_ee * Pz_ee + ss_Pxz_coef * Pz_ee * Px_ee + ss_Pxx_coef * Px_ee * Px_ee;
    OnDiagInfoType on_diag_info{
        // spin-first-site: g, spin-second-site: g
        {{gg, gg}, orbit_gg_gg - 0.25 * spinorbit_gg_gg - 0.25 * ss_coef},
        {{gg, ge}, orbit_ge_ge - 0.25 * spinorbit_ge_ge - 0.25 * ss_coef},
        {{ge, gg}, orbit_eg_eg - 0.25 * spinorbit_eg_eg - 0.25 * ss_coef},
        {{ge, ge}, orbit_ee_ee - 0.25 * spinorbit_ee_ee - 0.25 * ss_coef},
        // spin-first-site: g, spin-second-site: e
        {{gg, eg}, orbit_gg_gg + 0.25 * spinorbit_gg_gg + 0.25 * ss_coef},
        {{gg, ee}, orbit_ge_ge + 0.25 * spinorbit_ge_ge + 0.25 * ss_coef},
        {{ge, eg}, orbit_eg_eg + 0.25 * spinorbit_eg_eg + 0.25 * ss_coef},
        {{ge, ee}, orbit_ee_ee + 0.25 * spinorbit_ee_ee + 0.25 * ss_coef},
        // spin-first-site: e, spin-second-site: g
        {{eg, gg}, orbit_gg_gg + 0.25 * spinorbit_gg_gg + 0.25 * ss_coef},
        {{eg, ge}, orbit_ge_ge + 0.25 * spinorbit_ge_ge + 0.25 * ss_coef},
        {{ee, gg}, orbit_eg_eg + 0.25 * spinorbit_eg_eg + 0.25 * ss_coef},
        {{ee, ge}, orbit_ee_ee + 0.25 * spinorbit_ee_ee + 0.25 * ss_coef},
        // spin-first-site: e, spin-second-site: e
        {{eg, eg}, orbit_gg_gg - 0.25 * spinorbit_gg_gg - 0.25 * ss_coef},
        {{eg, ee}, orbit_ge_ge - 0.25 * spinorbit_ge_ge - 0.25 * ss_coef},
        {{ee, eg}, orbit_eg_eg - 0.25 * spinorbit_eg_eg - 0.25 * ss_coef},
        {{ee, ee}, orbit_ee_ee - 0.25 * spinorbit_ee_ee - 0.25 * ss_coef}};
    // Off-diag:
    const double orbit_gg_ge = +Pzz_coef * Pz_gg * Pz_ge + Pxz_coef * Px_gg * Pz_ge + Pxz_coef * Pz_gg * Px_ge + Pxx_coef * Px_gg * Px_ge;
    const double orbit_gg_eg = +Pzz_coef * Pz_ge * Pz_gg + Pxz_coef * Px_ge * Pz_gg + Pxz_coef * Pz_ge * Px_gg + Pxx_coef * Px_ge * Px_gg;
    const double orbit_gg_ee = +Pzz_coef * Pz_ge * Pz_ge + Pxz_coef * Px_ge * Pz_ge + Pxz_coef * Pz_ge * Px_ge + Pxx_coef * Px_ge * Px_ge;
    const double orbit_ge_eg = +Pzz_coef * Pz_ge * Pz_eg + Pxz_coef * Px_ge * Pz_eg + Pxz_coef * Pz_ge * Px_eg + Pxx_coef * Px_ge * Px_eg;
    const double orbit_ge_ee = +Pzz_coef * Pz_ge * Pz_ee + Pxz_coef * Px_ge * Pz_ee + Pxz_coef * Pz_ge * Px_ee + Pxx_coef * Px_ge * Px_ee;
    const double orbit_eg_ee = +Pzz_coef * Pz_ee * Pz_ge + Pxz_coef * Px_ee * Pz_ge + Pxz_coef * Pz_ee * Px_ge + Pxx_coef * Px_ee * Px_ge;
    const double spinorbit_gg_ge = +ss_Pzz_coef * Pz_gg * Pz_ge + ss_Pxz_coef * Px_gg * Pz_ge + ss_Pxz_coef * Pz_gg * Px_ge + ss_Pxx_coef * Px_gg * Px_ge;
    const double spinorbit_gg_eg = +ss_Pzz_coef * Pz_ge * Pz_gg + ss_Pxz_coef * Px_ge * Pz_gg + ss_Pxz_coef * Pz_ge * Px_gg + ss_Pxx_coef * Px_ge * Px_gg;
    const double spinorbit_gg_ee = +ss_Pzz_coef * Pz_ge * Pz_ge + ss_Pxz_coef * Px_ge * Pz_ge + ss_Pxz_coef * Pz_ge * Px_ge + ss_Pxx_coef * Px_ge * Px_ge;
    const double spinorbit_ge_eg = +ss_Pzz_coef * Pz_ge * Pz_eg + ss_Pxz_coef * Px_ge * Pz_eg + ss_Pxz_coef * Pz_ge * Px_eg + ss_Pxx_coef * Px_ge * Px_eg;
    const double spinorbit_ge_ee = +ss_Pzz_coef * Pz_ge * Pz_ee + ss_Pxz_coef * Px_ge * Pz_ee + ss_Pxz_coef * Pz_ge * Px_ee + ss_Pxx_coef * Px_ge * Px_ee;
    const double spinorbit_eg_ee = +ss_Pzz_coef * Pz_ee * Pz_ge + ss_Pxz_coef * Px_ee * Pz_ge + ss_Pxz_coef * Pz_ee * Px_ge + ss_Pxx_coef * Px_ee * Px_ge;
    OffDiagInfoType half_off_diag_info{
        // spin: not-changing, orbit: changing. // spin-first-site: g<->g, spin-second-site: g<->g, orbit: changing
        {{gg, gg}, {{gg, ge}, orbit_gg_ge - 0.25 * spinorbit_gg_ge}},
        {{gg, gg}, {{ge, gg}, orbit_gg_eg - 0.25 * spinorbit_gg_eg}},
        {{gg, gg}, {{ge, ge}, orbit_gg_ee - 0.25 * spinorbit_gg_ee}},
        {{gg, ge}, {{ge, gg}, orbit_ge_eg - 0.25 * spinorbit_ge_eg}},
        {{gg, ge}, {{ge, ge}, orbit_ge_ee - 0.25 * spinorbit_ge_ee}},
        {{ge, gg}, {{ge, ge}, orbit_eg_ee - 0.25 * spinorbit_eg_ee}},
        // spin: not-changing, orbit: changing. // spin-first-site: g<->g, spin-second-site: e<->e, orbit: changing
        {{gg, eg}, {{gg, ee}, orbit_gg_ge + 0.25 * spinorbit_gg_ge}},
        {{gg, eg}, {{ge, eg}, orbit_gg_eg + 0.25 * spinorbit_gg_eg}},
        {{gg, eg}, {{ge, ee}, orbit_gg_ee + 0.25 * spinorbit_gg_ee}},
        {{gg, ee}, {{ge, eg}, orbit_ge_eg + 0.25 * spinorbit_ge_eg}},
        {{gg, ee}, {{ge, ee}, orbit_ge_ee + 0.25 * spinorbit_ge_ee}},
        {{ge, eg}, {{ge, ee}, orbit_eg_ee + 0.25 * spinorbit_eg_ee}},
        // spin: not-changing, orbit: changing. // spin-first-site: e<->e, spin-second-site: g<->g, orbit: changing
        {{eg, gg}, {{eg, ge}, orbit_gg_ge + 0.25 * spinorbit_gg_ge}},
        {{eg, gg}, {{ee, gg}, orbit_gg_eg + 0.25 * spinorbit_gg_eg}},
        {{eg, gg}, {{ee, ge}, orbit_gg_ee + 0.25 * spinorbit_gg_ee}},
        {{eg, ge}, {{ee, gg}, orbit_ge_eg + 0.25 * spinorbit_ge_eg}},
        {{eg, ge}, {{ee, ge}, orbit_ge_ee + 0.25 * spinorbit_ge_ee}},
        {{ee, gg}, {{ee, ge}, orbit_eg_ee + 0.25 * spinorbit_eg_ee}},
        // spin: not-changing, orbit: changing. // spin-first-site: e<->e, spin-second-site: e<->e, orbit: changing
        {{eg, eg}, {{eg, ee}, orbit_gg_ge - 0.25 * spinorbit_gg_ge}},
        {{eg, eg}, {{ee, eg}, orbit_gg_eg - 0.25 * spinorbit_gg_eg}},
        {{eg, eg}, {{ee, ee}, orbit_gg_ee - 0.25 * spinorbit_gg_ee}},
        {{eg, ee}, {{ee, eg}, orbit_ge_eg - 0.25 * spinorbit_ge_eg}},
        {{eg, ee}, {{ee, ee}, orbit_ge_ee - 0.25 * spinorbit_ge_ee}},
        {{ee, eg}, {{ee, ee}, orbit_eg_ee - 0.25 * spinorbit_eg_ee}},
        // spin: changing, orbit: not-changing.
        {{gg, gg}, {{eg, eg}, 0.5 * ss_coef}},
        {{gg, ge}, {{eg, ee}, 0.5 * ss_coef}},
        {{ge, gg}, {{ee, eg}, 0.5 * ss_coef}},
        {{ge, ge}, {{ee, ee}, 0.5 * ss_coef}},
        // spin: changing, orbit: changing. // spin-first-site: g<->e, spin-second-site: g<->e.
        {{eg, eg}, {{gg, ge}, 0.5 * spinorbit_gg_ge}},
        {{eg, eg}, {{ge, gg}, 0.5 * spinorbit_gg_eg}},
        {{eg, eg}, {{ge, ge}, 0.5 * spinorbit_gg_ee}},
        {{eg, ee}, {{ge, gg}, 0.5 * spinorbit_ge_eg}},
        {{eg, ee}, {{ge, ge}, 0.5 * spinorbit_ge_ee}},
        {{ee, eg}, {{ge, ge}, 0.5 * spinorbit_eg_ee}},
        // spin: changing, orbit: changing. // spin-first-site: g<->e, spin-second-site: g<->e.
        {{gg, gg}, {{eg, ee}, 0.5 * spinorbit_gg_ge}},
        {{gg, gg}, {{ee, eg}, 0.5 * spinorbit_gg_eg}},
        {{gg, gg}, {{ee, ee}, 0.5 * spinorbit_gg_ee}},
        {{gg, ge}, {{ee, eg}, 0.5 * spinorbit_ge_eg}},
        {{gg, ge}, {{ee, ee}, 0.5 * spinorbit_ge_ee}},
        {{ge, gg}, {{ee, ee}, 0.5 * spinorbit_eg_ee}}};
    return chainkernel::OperatorKernel12<SoSiteStateTrait>{on_diag_info, half_off_diag_info};
}

chainkernel::OperatorKernel12<so_system::SoSiteStateTrait>
prepare_hamiltonian_kernel_12_fo(const HamiltonianParamsAfFo& params, double orbital_theta) {
    using namespace so_system;
    return prepare_hamiltonian_kernel_12_fo(
        params.get_ss_coef(),
        params.get_Pzz_coef(),
        params.get_Pxz_coef(),
        params.get_Pxx_coef(),
        params.get_ss_Pzz_coef(),
        params.get_ss_Pxz_coef(),
        params.get_ss_Pxx_coef(),
        orbital_theta);
}

chainkernel::OperatorKernel1<so_system::SoSiteStateTrait>
prepare_hamiltonian_kernel_1_fo(
    double s_coef,
    double tau_z_coef, double tau_minus_coef,
    double orbital_theta) {
    using namespace so_system;
    using OnDiagInfoType = std::map<chainkernel::StateKernel1<SoSiteStateTrait>, double>;
    using OffDiagInfoType = std::multimap<chainkernel::StateKernel1<SoSiteStateTrait>, chainkernel::CoupleInfoKernel1<SoSiteStateTrait>>;
    const double s1 = std::sin(orbital_theta), c1 = std::cos(orbital_theta);
    const double tau_z_gg = +c1;
    const double tau_z_ge = -s1;
    const double tau_z_ee = -c1;
    const double tau_minus_gg = -s1;
    const double tau_minus_ge = -c1;
    const double tau_minus_ee = +s1;
    const double orbit_g_g = tau_z_coef * tau_z_gg + tau_minus_coef * tau_minus_gg;
    const double orbit_e_e = tau_z_coef * tau_z_ee + tau_minus_coef * tau_minus_ee;
    const double spin_g_g = -0.5 * s_coef;
    const double spin_e_e = +0.5 * s_coef;
    OnDiagInfoType on_diag_info{
        {{gg}, spin_g_g + orbit_g_g},
        {{ge}, spin_g_g + orbit_e_e},
        {{eg}, spin_e_e + orbit_g_g},
        {{ee}, spin_e_e + orbit_e_e}};
    const double orbit_g_e = tau_z_coef * tau_z_ge + tau_minus_coef * tau_minus_ge;
    OffDiagInfoType half_off_diag_info{
        {{gg}, {{ge}, orbit_g_e}},
        {{eg}, {{ee}, orbit_g_e}}};
    return chainkernel::OperatorKernel1<SoSiteStateTrait>{on_diag_info, half_off_diag_info};
}

chainkernel::OperatorKernel1<so_system::SoSiteStateTrait>
prepare_hamiltonian_kernel_1_fo(const HamiltonianParamsAfFo& params, double orbital_theta) {
    using namespace so_system;
    return prepare_hamiltonian_kernel_1_fo(
        params.get_s_coef(),
        params.get_tau_z_coef(),
        params.get_tau_minus_coef(),
        orbital_theta);
}

}  // end of namespace so_hamiltonians
