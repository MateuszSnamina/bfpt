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

/*
 *
 * The below is for orbital hamiltonians defined for systems of `eg`-type orbitals
 * (where `eg` is a symbol of a irreducible representation
 *  of `Oh` point symmetry group).
 *
 * The relevant `eg` oribtal space consists of linear combinations of
 * the two basis oribtals: |x² - y²⟩ and |3z² - r²⟩ .
 *
 * Hamiltonians are often expressed in terms of
 * one-orbital projection operatos like these four:
 *
 * P^x = |x⟩⟨x|,         P^z = |x⟩⟨x|,
 * P^+ = |+⟩⟨+|,         P^- = |-⟩⟨-|,
 * P^Δ = (P^-) - (P^+)
 * where:
 *  |x⟩ and |z⟩ are shortcu notation of |x² - y²⟩ and |3z² - r²⟩,
 *  |+⟩ and |-⟩ denotes (|z⟩+|x⟩)/√2 and (|z⟩-|x⟩)/√2.
 *
 * The matrix representatrion of the projection operators
 * in the (|z⟩, |x⟩) basis is as follow:
 *
 *          |z⟩  |x⟩                          |z⟩  |x⟩
 *       ╭           ╮                     ╭           ╮
 * P^x = │  0    0   │  |z⟩          P^z = │ +1    0   │  |z⟩
 *       │  0   +1   │  |x⟩                │  0    0   │  |x⟩
 *       ╰           ╯                     ╰           ╯
 *
 *          |z⟩  |x⟩                          |z⟩  |x⟩
 *       ╭           ╮                     ╭           ╮
 * P^+ = │ +1/2 +1/2 │  |z⟩          P^- = │ +1/2 -1/2 │  |z⟩
 *       │ +1/2 +1/2 │  |x⟩                │ -1/2 +1/2 │  |x⟩
 *       ╰           ╯                     ╰           ╯
 *
 * In the program all calculations are performed in "optimal orbitals"
 * basis. It is defined by a single parameter θ_opt by the relation:
 *
 * ╭         ╮     ╭         ╮ ╭                                ╮
 * │ |∥⟩ |⟂⟩ │  =  │ |x⟩ |z⟩ │ │  +cos(θ_opt/2)   -sin(θ_opt/2) │
 * ╰         ╯     ╰         ╯ │  +sin(θ_opt/2)   +cos(θ_opt/2) │
 *                             ╰                                ╯
 *
 * Below there are expressions for the projection operators matrices
 * in "optimal orbitals" basis:
 *
 *          |∥⟩      |⟂⟩                      |∥⟩      |⟂⟩
 *       ╭                 ╮               ╭                 ╮
 * P^x = │ +(s2)²   +s2*c2 │ |∥⟩     P^z = │ +(c2)²   -s2*c2 │ |∥⟩
 *       │ +s2*c2   +(c2)² │ |⟂⟩           │ -s2*c2   +(s2)² │ |⟂⟩
 *       ╰                 ╯               ╰                 ╯
 * Where: c2 ≣ cos(θ_opt/2), s2 ≣ sin(θ_opt/2).
 *
 *                  |∥⟩  |⟂⟩                          |∥⟩   |⟂⟩
 *               ╭          ╮                      ╭           ╮
 * P^+ = ½ + ½ * │ +s1  +c1 │ |∥⟩    P^- = ½ + ½ * │ -s1   -c1 │ |∥⟩
 *               │ +c1  -s1 │ |⟂⟩                  │ -c1   +s1 │ |⟂⟩
 *               ╰          ╯                      ╰           ╯
 * Where: c1 ≣ cos(θ_opt), s1 ≣ sin(θ_opt).
 *
 *
 * The in the program the a general Hamiltonian for the orbitals chain is considered:
 * H = H_1 + H_12
 * H_1 = sum_i PΔ_coef * P^Δ(i)
 * H_12zz = sum_<ij> Pzz_coef * (P^z(i) P^z(j))
 * H_12xz = sum_<ij> Pxz_coef * (P^x(i) P^z(j) + P^z(i) P^x(j))
 * H_12xx = sum_<ij> Pxx_coef * (P^x(i) P^x(j))
 * H_12 = H_12zz + H_12xz + H_12xx
 */

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_kernel_12_af_fo(double Pzz_coef, double Pxz_coef, double Pxx_coef, double theta_opt) {
    using OnDiagInfoType = std::map<bfpt_common::StateKernel12<MonostarSiteState>, double>;
    using OffDiagInfoType = std::multimap<bfpt_common::StateKernel12<MonostarSiteState>, bfpt_common::CoupleInfoKernel12<MonostarSiteState>>;
    const double s2 = std::sin(theta_opt/2);
    const double c2 = std::cos(theta_opt/2);
    const double Px_gg = +s2 * s2;
    const double Px_ge = +s2 * c2;
    const double Px_eg = +s2 * c2;
    const double Px_ee = +c2 * c2;
    const double Pz_gg = +c2 * c2;
    const double Pz_ge = -s2 * c2;
    const double Pz_eg = -s2 * c2;
    const double Pz_ee = +s2 * s2;
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

bfpt_common::HamiltonianKernel1<MonostarSiteState>
prepare_hamiltonian_kernel_1_fo(double Pdelta_coef, double theta_opt) {
    using OnDiagInfoType = std::map<bfpt_common::StateKernel1<MonostarSiteState>, double>;
    using OffDiagInfoType = std::multimap<bfpt_common::StateKernel1<MonostarSiteState>, bfpt_common::CoupleInfoKernel1<MonostarSiteState>>;
    const double s1 = std::sin(theta_opt), c1 = std::cos(theta_opt);
    const double P_delta_gg = -s1;
    const double P_delta_ge = -c1;
    const double P_delta_ee = +s1;
    OnDiagInfoType on_diag_info {
        {{gs}, Pdelta_coef * P_delta_gg},
        {{es}, Pdelta_coef * P_delta_ee},
    };
    OffDiagInfoType half_off_diag_info {
        {{gs}, {{es}, Pdelta_coef * P_delta_ge}},
    };
    return bfpt_common::HamiltonianKernel1<MonostarSiteState>{on_diag_info, half_off_diag_info};
}

} // end of namespace model_monostar
