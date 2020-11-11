#pragma once

// #######################################################################
// ## HamiltonianParamsAfFm                                             ##
// #######################################################################

/*
 * In case od AF:
 * H₁ = - B * Σ_i^onLatticeA Sᶻ(i) + B * Σ_i^onLatticeB Sᶻ(i)
 * H₁₂ = + J_classical Σ_<ij> (Sᶻ(i) Sᶻ(j)) + J_quantum * 0.5 * Σ_<ij> [S⁺(i) S⁻(j) + S⁻(i) S⁺(j)]
 * H = H₁ + H₁₂ + free_coef * n_sites
 *
 * onLatticeA: |gs⟩ = |↑⟩, |es⟩ = |↓⟩
 * onLatticeB: |gs⟩ = |↓⟩, |es⟩ = |↑⟩
 *
 *
 * In case od FM:
 * H₁ = - B * Σ_i Sᶻ(i)
 * H₁₂ = - J_classical Σ_<ij> (Sᶻ(i) Sᶻ(j)) - J_quantum * 0.5 * Σ_<ij> [S⁺(i) S⁻(j) + S⁻(i) S⁺(j)]
 * H = H₁ + H₁₂ + free_coef * n_sites
 *
 * on each site: |gs⟩ = |↑⟩, |es> = |↓⟩
 */

namespace monostar_hamiltonians {

class HamiltonianParamsAfFm {
   public:
    class Builder {
       public:
        Builder set_J_classical(double);
        Builder set_J_quantum(double);
        Builder set_B(double);
        Builder set_free(double);
        HamiltonianParamsAfFm build() const;
       private:
        double _J_classicall = 1.0;
        double _J_quantum = 1.0;
        double _B = 0.0;
        double _free = 0.0;
    };
    // Construction destruction:
    friend HamiltonianParamsAfFm Builder::build() const;
    // Getters:
    double get_J_classical() const;
    double get_J_quantum() const;
    double get_B() const;
    double get_free() const;
   private:
    HamiltonianParamsAfFm(double J_classicall, double J_quantum, double B, double free);
    double _J_classicall;
    double _J_quantum;
    double _B;
    double _free;
};

}  // end of namespace monostar_hamiltonians
