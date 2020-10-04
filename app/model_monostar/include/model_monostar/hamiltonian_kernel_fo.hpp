#pragma once

#include<model_monostar/monostar_site_state.hpp>
#include<model_monostar/hamiltonian_params_fo.hpp>

#include<bfpt_common/operator_kernel.hpp>

#include<map>

// #######################################################################
// ## prepare_hamiltonian_kernel_{1,12}_fo                              ##
// #######################################################################

/*
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
 * Pᶻ = |z⟩⟨z|,         Pˣ = |x⟩⟨x|,
 * P⁺ = |+⟩⟨+|,         P⁻ = |-⟩⟨-|,
 *
 * where:
 *  |x⟩ and |z⟩ are shortcu notation of |x² - y²⟩ and |3z² - r²⟩,
 *  |+⟩ and |-⟩ denotes (|z⟩+|x⟩)/√2 and (|z⟩-|x⟩)/√2.
 *
 * The matrix representatrion of the projection operators
 * in the (|z⟩, |x⟩) basis is as follow:
 *
 *        |z⟩ |x⟩                       |z⟩ |x⟩
 *      ╭        ╮                    ╭        ╮
 * Pᶻ = │ +1   0 │  |z⟩          Pˣ = │  0   0 │  |z⟩
 *      │  0   0 │  |x⟩               │  0  +1 │  |x⟩
 *      ╰        ╯                    ╰        ╯
 *
 *        |z⟩ |x⟩                       |z⟩ |x⟩
 *      ╭        ╮                    ╭        ╮
 * P⁺ = │ +½  +½ │  |z⟩          P⁻ = │ +½  -½ │  |z⟩
 *      │ +½  +½ │  |x⟩               │ -½  +½ │  |x⟩
 *      ╰        ╯                    ╰        ╯
 *
 * In the program all calculations are performed in "rotated orbitals"
 * basis. It is defined by a single parameter θ (called "orbital_theta")
 * by the relation:
 *
 * ╭         ╮     ╭         ╮ ╭                        ╮
 * │ |∥⟩ |⟂⟩ │  =  │ |x⟩ |z⟩ │ │  +cos(θ/2)   -sin(θ/2) │
 * ╰         ╯     ╰         ╯ │  +sin(θ/2)   +cos(θ/2) │
 *                             ╰                        ╯
 *
 * Below there are expressions for the projection operators matrices
 * in "rotated orbitals" basis:
 *
 *                 |∥⟩  |⟂⟩            |∥⟩      |⟂⟩
 *              ╭          ╮       ╭                 ╮
 * Pᶻ = ½ + ½ * │ +c1  -s1 │ |∥⟩ = │ +(c2)²   -s2*c2 │ |∥⟩
 *              │ -s1  -c1 │ |⟂⟩   │ -s2*c2   +(s2)² │ |⟂⟩
 *              ╰          ╯       ╰                 ╯
 *
 *                 |∥⟩  |⟂⟩            |∥⟩      |⟂⟩
 *              ╭          ╮       ╭                 ╮
 * Pˣ = ½ + ½ * │ -c1  +s1 │ |∥⟩ = │ +(s2)²   +s2*c2 │ |∥⟩
 *              │ +s1  +c1 │ |⟂⟩   │ +s2*c2   +(c2)² │ |⟂⟩
 *              ╰          ╯       ╰                 ╯
 *
 *                 |∥⟩  |⟂⟩               |∥⟩         |⟂⟩
 *              ╭          ╮           ╭                         ╮
 * P⁺ = ½ + ½ * │ +s1  +c1 │ |∥⟩ = ½ * │ 1+2*s2*c2   (c2)²-(s2)² │ |∥⟩
 *              │ +c1  -s1 │ |⟂⟩       │ (c2)²-(s2)² 1-2*s2*c2   │ |⟂⟩
 *              ╰          ╯           ╰                         ╯
 *
 *                |∥⟩  |⟂⟩               |∥⟩         |⟂⟩
 *              ╭          ╮           ╭                         ╮
 * P⁻ = ½ + ½ * │ -s1  -c1 │ |∥⟩ = ½ * │ 1-2*s2*c2   (s2)²-(c2)² │ |∥⟩
 *              │ -c1  +s1 │ |⟂⟩       │ (s2)²-(c2)² 1+2*s2*c2   │ |⟂⟩
 *              ╰          ╯           ╰                         ╯
 *
 * Where: c2 ≣ cos(θ/2), s2 ≣ sin(θ/2),
 *        c1 ≣ cos(θ),   s1 ≣ sin(θ).
 *
 * We also define τ operators:
 *
 * τᶻ = |z⟩⟨z| - |x⟩⟨x|
 * τ⁻ = |-⟩⟨-| - |+⟩⟨+|
 *
 * The in the program the a general Hamiltonian for the orbitals chain is considered:
 * H₁ᶻ = sum_i tau_z_coef * τᶻ(i)
 * H₁⁻ = sum_i tau_minus_coef * τ⁻(i)
 * H₁₂ᶻᶻ = sum_<ij> Pzz_coef * (Pᶻ(i) Pᶻ(j))
 * H₁₂ˣᶻ = sum_<ij> Pxz_coef * (Pˣ(i) Pᶻ(j) + Pᶻ(i) Pˣ(j))
 * H₁₂ˣˣ = sum_<ij> Pxx_coef * (Pˣ(i) Pˣ(j))
 * H₁ = H₁⁻ + H₁ᶻ
 * H₁₂ = H₁₂ᶻᶻ + H₁₂ˣᶻ + H₁₂ˣˣ
 * H = H₁ + H₁₂
 */

namespace model_monostar {

bfpt_common::OperatorKernel12<MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_fo(double Pzz_coef, double Pxz_coef, double Pxx_coef, double orbital_theta);

bfpt_common::OperatorKernel12<MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_fo(const HamiltonianParamsFo&, double orbital_theta);

bfpt_common::OperatorKernel1<MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1_fo(double tau_z_coef, double tau_minus_coef, double orbital_theta);

bfpt_common::OperatorKernel1<MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1_fo(const HamiltonianParamsFo&, double orbital_theta);

} // end of namespace model_monostar
