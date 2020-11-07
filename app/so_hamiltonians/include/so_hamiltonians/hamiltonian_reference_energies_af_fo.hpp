#pragma once

#include <monostar_hamiltonians/hamiltonian_reference_energies.hpp>
#include <so_hamiltonians/hamiltonian_params_af_fo.hpp>

// #######################################################################
// ## HamiltonianReferenceEnergiesFo                                    ##
// #######################################################################

namespace so_hamiltonians {

class HamiltonianReferenceEnergiesFo final : public monostar_hamiltonians::HamiltonianReferenceEnergies {
   public:
    HamiltonianReferenceEnergiesFo(unsigned n_sites, const HamiltonianParamsAfFo& params, double orbital_theta);
    std::optional<double> get_gs_energy() const override;
    std::optional<double> get_es_exciation_enery(unsigned n_k) const override;

   private:
    const HamiltonianParamsAfFo& _params;
    const double _orbital_theta;
};

}  // end of namespace so_hamiltonians
