#include <so_hamiltonians/hamiltonian_kernel_af_fo.hpp>

#include <monostar_hamiltonians/hamiltonian_params_fo_site_matrices.hpp>

// #######################################################################
// ## prepare_hamiltonian_kernel_{1,12}_af_fo                           ##
// #######################################################################

namespace so_hamiltonians {

chainkernel::OperatorKernel12<so_system::SoSiteStateTrait>
prepare_hamiltonian_kernel_12_af_fo(
    double ss_coef,
    double Pzz_coef, double Pxz_coef, double Pxx_coef,
    double ss_Pzz_coef, double ss_Pxz_coef, double ss_Pxx_coef,
    double orbital_theta) {
    using namespace so_system;
    using OnDiagInfoType = std::map<chainkernel::StateKernel12<SoSiteStateTrait>, double>;
    using OffDiagInfoType = std::multimap<chainkernel::StateKernel12<SoSiteStateTrait>, chainkernel::CoupleInfoKernel12<SoSiteStateTrait>>;
    const double Pz_gg = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta)(0, 0));
    const double Pz_ge = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta)(0, 1));
    const double Pz_eg = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta)(1, 0));
    const double Pz_ee = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_z_in_ge_basis(orbital_theta)(1, 1));
    const double Px_gg = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta)(0, 0));
    const double Px_ge = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta)(0, 1));
    const double Px_eg = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta)(1, 0));
    const double Px_ee = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_P_x_in_ge_basis(orbital_theta)(1, 1));
    // On-diag:
    const double pureorbit_gg_gg = +Pzz_coef * Pz_gg * Pz_gg + Pxz_coef * Px_gg * Pz_gg + Pxz_coef * Pz_gg * Px_gg + Pxx_coef * Px_gg * Px_gg;
    const double pureorbit_ge_ge = +Pzz_coef * Pz_gg * Pz_ee + Pxz_coef * Px_gg * Pz_ee + Pxz_coef * Pz_gg * Px_ee + Pxx_coef * Px_gg * Px_ee;
    const double pureorbit_eg_eg = +Pzz_coef * Pz_ee * Pz_gg + Pxz_coef * Px_ee * Pz_gg + Pxz_coef * Pz_ee * Px_gg + Pxx_coef * Px_ee * Px_gg;
    const double pureorbit_ee_ee = +Pzz_coef * Pz_ee * Pz_ee + Pxz_coef * Px_ee * Pz_ee + Pxz_coef * Pz_ee * Px_ee + Pxx_coef * Px_ee * Px_ee;
    const double spinorbit_gg_gg = +ss_Pzz_coef * Pz_gg * Pz_gg + ss_Pxz_coef * Px_gg * Pz_gg + ss_Pxz_coef * Pz_gg * Px_gg + ss_Pxx_coef * Px_gg * Px_gg;
    const double spinorbit_ge_ge = +ss_Pzz_coef * Pz_gg * Pz_ee + ss_Pxz_coef * Px_gg * Pz_ee + ss_Pxz_coef * Pz_gg * Px_ee + ss_Pxx_coef * Px_gg * Px_ee;
    const double spinorbit_eg_eg = +ss_Pzz_coef * Pz_ee * Pz_gg + ss_Pxz_coef * Px_ee * Pz_gg + ss_Pxz_coef * Pz_ee * Px_gg + ss_Pxx_coef * Px_ee * Px_gg;
    const double spinorbit_ee_ee = +ss_Pzz_coef * Pz_ee * Pz_ee + ss_Pxz_coef * Px_ee * Pz_ee + ss_Pxz_coef * Pz_ee * Px_ee + ss_Pxx_coef * Px_ee * Px_ee;
    OnDiagInfoType on_diag_info{
        // spin-first-site: g, spin-second-site: g
        {{gg, gg}, pureorbit_gg_gg - 0.25 * spinorbit_gg_gg - 0.25 * ss_coef},
        {{gg, ge}, pureorbit_ge_ge - 0.25 * spinorbit_ge_ge - 0.25 * ss_coef},
        {{ge, gg}, pureorbit_eg_eg - 0.25 * spinorbit_eg_eg - 0.25 * ss_coef},
        {{ge, ge}, pureorbit_ee_ee - 0.25 * spinorbit_ee_ee - 0.25 * ss_coef},
        // spin-first-site: g, spin-second-site: e
        {{gg, eg}, pureorbit_gg_gg + 0.25 * spinorbit_gg_gg + 0.25 * ss_coef},
        {{gg, ee}, pureorbit_ge_ge + 0.25 * spinorbit_ge_ge + 0.25 * ss_coef},
        {{ge, eg}, pureorbit_eg_eg + 0.25 * spinorbit_eg_eg + 0.25 * ss_coef},
        {{ge, ee}, pureorbit_ee_ee + 0.25 * spinorbit_ee_ee + 0.25 * ss_coef},
        // spin-first-site: e, spin-second-site: g
        {{eg, gg}, pureorbit_gg_gg + 0.25 * spinorbit_gg_gg + 0.25 * ss_coef},
        {{eg, ge}, pureorbit_ge_ge + 0.25 * spinorbit_ge_ge + 0.25 * ss_coef},
        {{ee, gg}, pureorbit_eg_eg + 0.25 * spinorbit_eg_eg + 0.25 * ss_coef},
        {{ee, ge}, pureorbit_ee_ee + 0.25 * spinorbit_ee_ee + 0.25 * ss_coef},
        // spin-first-site: e, spin-second-site: e
        {{eg, eg}, pureorbit_gg_gg - 0.25 * spinorbit_gg_gg - 0.25 * ss_coef},
        {{eg, ee}, pureorbit_ge_ge - 0.25 * spinorbit_ge_ge - 0.25 * ss_coef},
        {{ee, eg}, pureorbit_eg_eg - 0.25 * spinorbit_eg_eg - 0.25 * ss_coef},
        {{ee, ee}, pureorbit_ee_ee - 0.25 * spinorbit_ee_ee - 0.25 * ss_coef}};
    // Off-diag:
    const double pureorbit_gg_ge = +Pzz_coef * Pz_gg * Pz_ge + Pxz_coef * Px_gg * Pz_ge + Pxz_coef * Pz_gg * Px_ge + Pxx_coef * Px_gg * Px_ge;
    const double pureorbit_gg_eg = +Pzz_coef * Pz_ge * Pz_gg + Pxz_coef * Px_ge * Pz_gg + Pxz_coef * Pz_ge * Px_gg + Pxx_coef * Px_ge * Px_gg;
    const double pureorbit_gg_ee = +Pzz_coef * Pz_ge * Pz_ge + Pxz_coef * Px_ge * Pz_ge + Pxz_coef * Pz_ge * Px_ge + Pxx_coef * Px_ge * Px_ge;
    const double pureorbit_ge_eg = +Pzz_coef * Pz_ge * Pz_eg + Pxz_coef * Px_ge * Pz_eg + Pxz_coef * Pz_ge * Px_eg + Pxx_coef * Px_ge * Px_eg;
    const double pureorbit_ge_ee = +Pzz_coef * Pz_ge * Pz_ee + Pxz_coef * Px_ge * Pz_ee + Pxz_coef * Pz_ge * Px_ee + Pxx_coef * Px_ge * Px_ee;
    const double pureorbit_eg_ee = +Pzz_coef * Pz_ee * Pz_ge + Pxz_coef * Px_ee * Pz_ge + Pxz_coef * Pz_ee * Px_ge + Pxx_coef * Px_ee * Px_ge;
    const double spinorbit_gg_ge = +ss_Pzz_coef * Pz_gg * Pz_ge + ss_Pxz_coef * Px_gg * Pz_ge + ss_Pxz_coef * Pz_gg * Px_ge + ss_Pxx_coef * Px_gg * Px_ge;
    const double spinorbit_gg_eg = +ss_Pzz_coef * Pz_ge * Pz_gg + ss_Pxz_coef * Px_ge * Pz_gg + ss_Pxz_coef * Pz_ge * Px_gg + ss_Pxx_coef * Px_ge * Px_gg;
    const double spinorbit_gg_ee = +ss_Pzz_coef * Pz_ge * Pz_ge + ss_Pxz_coef * Px_ge * Pz_ge + ss_Pxz_coef * Pz_ge * Px_ge + ss_Pxx_coef * Px_ge * Px_ge;
    const double spinorbit_ge_eg = +ss_Pzz_coef * Pz_ge * Pz_eg + ss_Pxz_coef * Px_ge * Pz_eg + ss_Pxz_coef * Pz_ge * Px_eg + ss_Pxx_coef * Px_ge * Px_eg;
    const double spinorbit_ge_ee = +ss_Pzz_coef * Pz_ge * Pz_ee + ss_Pxz_coef * Px_ge * Pz_ee + ss_Pxz_coef * Pz_ge * Px_ee + ss_Pxx_coef * Px_ge * Px_ee;
    const double spinorbit_eg_ee = +ss_Pzz_coef * Pz_ee * Pz_ge + ss_Pxz_coef * Px_ee * Pz_ge + ss_Pxz_coef * Pz_ee * Px_ge + ss_Pxx_coef * Px_ee * Px_ge;
    std::cout << "spinorbit_gg_ge: " << spinorbit_gg_ge << std::endl;  //TODO remove
    std::cout << "spinorbit_gg_eg: " << spinorbit_gg_eg << std::endl;  //TODO remove
    std::cout << "spinorbit_gg_ee: " << spinorbit_gg_ee << std::endl;  //TODO remove
    std::cout << "spinorbit_ge_eg: " << spinorbit_ge_eg << std::endl;  //TODO remove
    std::cout << "spinorbit_ge_ee: " << spinorbit_ge_ee << std::endl;  //TODO remove
    std::cout << "spinorbit_eg_ee: " << spinorbit_eg_ee << std::endl;  //TODO remove
    OffDiagInfoType half_off_diag_info{
        // purespin. spins: gg<->ee, orbitals: fixed.
        {{gg, gg}, {{eg, eg}, 0.5 * ss_coef + 0.5 * spinorbit_gg_gg}},  // with orbitals: gg
        {{gg, ge}, {{eg, ee}, 0.5 * ss_coef + 0.5 * spinorbit_ge_ge}},  // with orbitals: ge
        {{ge, gg}, {{ee, eg}, 0.5 * ss_coef + 0.5 * spinorbit_eg_eg}},  // with orbitals: eg
        {{ge, ge}, {{ee, ee}, 0.5 * ss_coef + 0.5 * spinorbit_ee_ee}},  // with orbitals: ee
        // pureorbit + spinorbit. orbitals: gg<->ge, spins: fixed.
        {{gg, gg}, {{gg, ge}, pureorbit_gg_ge - 0.25 * spinorbit_gg_ge}},  // with spins: gg
        {{gg, eg}, {{gg, ee}, pureorbit_gg_ge + 0.25 * spinorbit_gg_ge}},  // with spins: ge
        {{eg, gg}, {{eg, ge}, pureorbit_gg_ge + 0.25 * spinorbit_gg_ge}},  // with spins: eg
        {{eg, eg}, {{eg, ee}, pureorbit_gg_ge - 0.25 * spinorbit_gg_ge}},  // with spins: ee
        // pureorbit + spinorbit. orbitals: gg<->eg, spins: fixed.
        {{gg, gg}, {{ge, gg}, pureorbit_gg_eg - 0.25 * spinorbit_gg_eg}},  // with spins: gg
        {{gg, eg}, {{ge, eg}, pureorbit_gg_eg + 0.25 * spinorbit_gg_eg}},  // with spins: ge
        {{eg, gg}, {{ee, gg}, pureorbit_gg_eg + 0.25 * spinorbit_gg_eg}},  // with spins: eg
        {{eg, eg}, {{ee, eg}, pureorbit_gg_eg - 0.25 * spinorbit_gg_eg}},  // with spins: ee
        // pureorbit + spinorbit. orbitals: gg<->ee, spins: fixed.
        {{gg, gg}, {{ge, ge}, pureorbit_gg_ee - 0.25 * spinorbit_gg_ee}},  // with spins: gg
        {{gg, eg}, {{ge, ee}, pureorbit_gg_ee + 0.25 * spinorbit_gg_ee}},  // with spins: ge
        {{eg, gg}, {{ee, ge}, pureorbit_gg_ee + 0.25 * spinorbit_gg_ee}},  // with spins: eg
        {{eg, eg}, {{ee, ee}, pureorbit_gg_ee - 0.25 * spinorbit_gg_ee}},  // with spins: ee
        // pureorbit + spinorbit. orbitals: ge<->eg, spins: fixed.
        {{gg, ge}, {{ge, gg}, pureorbit_ge_eg - 0.25 * spinorbit_ge_eg}},  // with spins: gg
        {{gg, ee}, {{ge, eg}, pureorbit_ge_eg + 0.25 * spinorbit_ge_eg}},  // with spins: ge
        {{eg, ge}, {{ee, gg}, pureorbit_ge_eg + 0.25 * spinorbit_ge_eg}},  // with spins: eg
        {{eg, ee}, {{ee, eg}, pureorbit_ge_eg - 0.25 * spinorbit_ge_eg}},  // with spins: ee
        // pureorbit + spinorbit. orbitals: ge<->ee, spins: fixed.
        {{gg, ge}, {{ge, ge}, pureorbit_ge_ee - 0.25 * spinorbit_ge_ee}},  // with spins: gg
        {{gg, ee}, {{ge, ee}, pureorbit_ge_ee + 0.25 * spinorbit_ge_ee}},  // with spins: ge
        {{eg, ge}, {{ee, ge}, pureorbit_ge_ee + 0.25 * spinorbit_ge_ee}},  // with spins: eg
        {{eg, ee}, {{ee, ee}, pureorbit_ge_ee - 0.25 * spinorbit_ge_ee}},  // with spins: ee
        // pureorbit + spinorbit. orbitals: eg<->ee, spins: fixed.
        {{ge, gg}, {{ge, ge}, pureorbit_eg_ee - 0.25 * spinorbit_eg_ee}},  // with spins: gg
        {{ge, eg}, {{ge, ee}, pureorbit_eg_ee + 0.25 * spinorbit_eg_ee}},  // with spins: ge
        {{ee, gg}, {{ee, ge}, pureorbit_eg_ee + 0.25 * spinorbit_eg_ee}},  // with spins: eg
        {{ee, eg}, {{ee, ee}, pureorbit_eg_ee - 0.25 * spinorbit_eg_ee}},  // with spins: ee
        // spinorbit. spins: ee<->gg, orbit: changing.
        {{eg, eg}, {{gg, ge}, 0.5 * spinorbit_gg_ge}},
        {{eg, eg}, {{ge, gg}, 0.5 * spinorbit_gg_eg}},
        {{eg, eg}, {{ge, ge}, 0.5 * spinorbit_gg_ee}},
        {{eg, ee}, {{ge, gg}, 0.5 * spinorbit_ge_eg}},
        {{eg, ee}, {{ge, ge}, 0.5 * spinorbit_ge_ee}},
        {{ee, eg}, {{ge, ge}, 0.5 * spinorbit_eg_ee}},
        // spinorbit. spins: gg<->ee, orbit: changing.
        {{gg, gg}, {{eg, ee}, 0.5 * spinorbit_gg_ge}},
        {{gg, gg}, {{ee, eg}, 0.5 * spinorbit_gg_eg}},
        {{gg, gg}, {{ee, ee}, 0.5 * spinorbit_gg_ee}},
        {{gg, ge}, {{ee, eg}, 0.5 * spinorbit_ge_eg}},
        {{gg, ge}, {{ee, ee}, 0.5 * spinorbit_ge_ee}},
        {{ge, gg}, {{ee, ee}, 0.5 * spinorbit_eg_ee}}};
    return chainkernel::OperatorKernel12<SoSiteStateTrait>{on_diag_info, half_off_diag_info};
}

chainkernel::OperatorKernel12<so_system::SoSiteStateTrait>
prepare_hamiltonian_kernel_12_af_fo(const HamiltonianParamsAfFo& params, double orbital_theta) {
    using namespace so_system;
    return prepare_hamiltonian_kernel_12_af_fo(
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
prepare_hamiltonian_kernel_1_af_fo(
    double s_coef,
    double tau_z_coef, double tau_minus_coef,
    double orbital_theta) {
    using namespace so_system;
    using OnDiagInfoType = std::map<chainkernel::StateKernel1<SoSiteStateTrait>, double>;
    using OffDiagInfoType = std::multimap<chainkernel::StateKernel1<SoSiteStateTrait>, chainkernel::CoupleInfoKernel1<SoSiteStateTrait>>;
    const double tau_z_gg = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_z_in_ge_basis(orbital_theta)(0, 0));
    const double tau_z_ge = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_z_in_ge_basis(orbital_theta)(0, 1));
    const double tau_z_ee = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_z_in_ge_basis(orbital_theta)(1, 1));
    const double tau_minus_gg = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_minus_in_ge_basis(orbital_theta)(0, 0));
    const double tau_minus_ge = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_minus_in_ge_basis(orbital_theta)(0, 1));
    const double tau_minus_ee = std::real(monostar_hamiltonians::OneSiteOrbitalMatrices::get_tau_minus_in_ge_basis(orbital_theta)(1, 1));
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
prepare_hamiltonian_kernel_1_af_fo(const HamiltonianParamsAfFo& params, double orbital_theta) {
    using namespace so_system;
    return prepare_hamiltonian_kernel_1_af_fo(
        params.get_s_coef(),
        params.get_tau_z_coef(),
        params.get_tau_minus_coef(),
        orbital_theta);
}

}  // end of namespace so_hamiltonians
