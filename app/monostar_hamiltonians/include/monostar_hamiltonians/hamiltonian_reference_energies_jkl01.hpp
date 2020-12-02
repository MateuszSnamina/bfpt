#pragma once

#include <monostar_hamiltonians/hamiltonian_reference_energies.hpp>
#include <monostar_hamiltonians/hamiltonian_params_jkl01.hpp>

// #######################################################################
// ## HamiltonianReferenceEnergiesJkl01                                 ##
// #######################################################################

namespace monostar_hamiltonians {

class HamiltonianReferenceEnergiesJkl01 final : public HamiltonianReferenceEnergies {
   public:
    HamiltonianReferenceEnergiesJkl01(unsigned n_sites, const HamiltonianParamsJkl01& params);
    std::optional<double> get_gs_energy() const override;
    std::optional<double> get_es_exciation_enery(unsigned n_k) const override;

   private:
    const HamiltonianParamsJkl01& _params;
};

}  // end of namespace monostar_hamiltonians
