#pragma once

#include<model_monostar/hamiltonian_reference_energies.hpp>
#include<model_monostar/hamiltonian_params_fo.hpp>

// #######################################################################
// ## HamiltonianReferenceEnergiesFo                                    ##
// #######################################################################

namespace model_monostar {

class HamiltonianReferenceEnergiesFo final : public HamiltonianReferenceEnergies {
public:
    HamiltonianReferenceEnergiesFo(unsigned n_sites, const HamiltonianParamsFo& params, double orbital_theta);
    std::optional<double> get_gs_energy() const override;
    std::optional<double> get_es_exciation_enery(unsigned n_k) const override;
private:
    const HamiltonianParamsFo& _params;
    const double _orbital_theta;
};

} // end of namespace model_monostar
