#pragma once

#include<monostar_hamiltonians/hamiltonian_reference_energies.hpp>
#include<monostar_hamiltonians/hamiltonian_params_fo.hpp>

// #######################################################################
// ## HamiltonianReferenceEnergiesFo                                    ##
// #######################################################################

namespace monostar_hamiltonians {

class HamiltonianReferenceEnergiesFo final : public HamiltonianReferenceEnergies {
public:
    HamiltonianReferenceEnergiesFo(unsigned n_sites, const HamiltonianParamsFo& params, double orbital_theta);
    std::optional<double> get_gs_energy() const override;
    std::optional<double> get_es_exciation_enery(unsigned n_k) const override;
private:
    const HamiltonianParamsFo& _params;
    const double _orbital_theta;
};

} // end of namespace monostar_hamiltonians
