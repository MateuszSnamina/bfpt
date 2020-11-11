#pragma once

#include <monostar_hamiltonians/hamiltonian_reference_energies.hpp>
#include <monostar_hamiltonians/hamiltonian_params_af_fm.hpp>

#include <armadillo>

#include <optional>

// #######################################################################
// ## HamiltonianReferenceEnergiesFm                                    ##
// #######################################################################

namespace monostar_hamiltonians {

class HamiltonianReferenceEnergiesFm final : public HamiltonianReferenceEnergies {
   public:
    HamiltonianReferenceEnergiesFm(unsigned n_sites, double J_classical, double J_quantum, double B, double free);
    HamiltonianReferenceEnergiesFm(unsigned n_sites, const HamiltonianParamsAfFm& params);
    std::optional<double> get_gs_energy() const override;
    std::optional<double> get_es_exciation_enery(unsigned n_k) const override;

   private:
    const double _J_classical;
    const double _J_quantum;
    const double _B;
    const double _free;
};

}  // end of namespace monostar_hamiltonians

// #######################################################################
// ## HamiltonianReferenceEnergiesAf                                    ##
// #######################################################################

namespace monostar_hamiltonians {

class HamiltonianReferenceEnergiesAf final : public HamiltonianReferenceEnergies {
   public:
    HamiltonianReferenceEnergiesAf(unsigned n_sites, double J, double free);
    HamiltonianReferenceEnergiesAf(unsigned n_sites, const HamiltonianParamsAfFm& params);
    std::optional<double> get_gs_energy() const override;
    std::optional<double> get_es_exciation_enery(unsigned n_k) const override;

   private:
    const double _J;
    const double _free;
    const bool _is_applicable;
};

}  // end of namespace monostar_hamiltonians
