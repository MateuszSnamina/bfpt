#pragma once

#include <monostar_hamiltonians/hamiltonian_params_fo_helpers.hpp>

#include <set>
#include <string>

/*
 * The below is for orbital hamiltonians defined for systems of `e_g`-type orbitals
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
 *  |x⟩ and |z⟩ are shortcut notation of |x² - y²⟩ and |3z² - r²⟩,
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
 * ╭         ╮     ╭         ╮              ╭                        ╮
 * │ |∥⟩ |⟂⟩ │  =  │ |x⟩ |z⟩ │ β, where β = │  +cos(θ/2)   -sin(θ/2) │
 * ╰         ╯     ╰         ╯              │  +sin(θ/2)   +cos(θ/2) │
 *                                          ╰                        ╯
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
 * H₁ᶻ = Σ_i tau_z_coef * τᶻ(i)
 * H₁⁻ = Σ_i tau_minus_coef * τ⁻(i)
 * H₁₂ᶻᶻ = Σ_<ij> Pzz_coef * (Pᶻ(i) Pᶻ(j))
 * H₁₂ˣᶻ = Σ_<ij> Pxz_coef * (Pˣ(i) Pᶻ(j) + Pᶻ(i) Pˣ(j))
 * H₁₂ˣˣ = Σ_<ij> Pxx_coef * (Pˣ(i) Pˣ(j))
 * H₁ = H₁ᶻ + H₁⁻
 * H₁₂ = H₁₂ᶻᶻ + H₁₂ˣᶻ + H₁₂ˣˣ
 * H = H₁ + H₁₂ + free_coef * n_sites
 *
 */

// #######################################################################
// ## HamiltonianParamsFo                                               ##
// #######################################################################

namespace monostar_hamiltonians {

class HamiltonianParamsFo {
   public:
    // Builder:
    class Builder {
       public:
        Builder set_tau_z_coef(double);
        Builder set_tau_minus_coef(double);
        Builder set_Pzz_coef(double);
        Builder set_Pxz_coef(double);
        Builder set_Pxx_coef(double);
        Builder set_free_coef(double);
        HamiltonianParamsFo build() const;

       private:
        double _tau_z_coef = 0.0;
        double _tau_munis_coef = 0.0;
        double _Pzz_coef = 0.0;
        double _Pxz_coef = 0.0;
        double _Pxx_coef = 0.0;
        double _free_coef = 0.0;
    };
    // Construction/destruction:
    friend HamiltonianParamsFo Builder::build() const;
    // Getters:
    double get_tau_z_coef() const;
    double get_tau_minus_coef() const;
    double get_Pzz_coef() const;
    double get_Pxz_coef() const;
    double get_Pxx_coef() const;
    double get_free_coef() const;
    double get_site_energy(double theta) const;
    double get_site_energy_derivative(double theta) const;
    double get_site_energy_derivative2(double theta) const;
    double get_site_energy_derivative3(double theta) const;
    double get_site_energy_derivative4(double theta) const;
    std::set<double> get_theta_opt() const;
    std::set<double> get_theta_opt_numerical() const;
    utility::Result<std::set<double>, NoKnownAnalyticalSolutionError> get_theta_opt_analytical() const;
    std::string string_repr_in_orbital_operators() const;
    std::string string_repr_in_trigonometric_functions() const;
   private:
    HamiltonianParamsFo(double tau_z_coef, double tau_minus_coef, double Pzz_coef, double Pxz_coef, double Pxx_coef, double free_coef);
    double _tau_z_coef;
    double _tau_minus_coef;
    double _Pzz_coef;
    double _Pxz_coef;
    double _Pxx_coef;
    double _free_coef = 0.0;
};

}  // end of namespace monostar_hamiltonians
