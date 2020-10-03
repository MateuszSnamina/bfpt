#include<model_monostar/hamiltonian_reference_energies_fo.hpp>

#include<cassert>

// #######################################################################
// ## HamiltonianReferenceEnergiesFo                                    ##
// #######################################################################

namespace model_monostar {

HamiltonianReferenceEnergiesFo::HamiltonianReferenceEnergiesFo(unsigned n_sites, const HamiltonianParamsFo& params, double orbital_theta) :
    HamiltonianReferenceEnergies(n_sites),
    _params(params),
  _orbital_theta(orbital_theta) {
}

std::optional<double> HamiltonianReferenceEnergiesFo::get_gs_energy() const {
    return _n_sites * _params.get_site_energy(_orbital_theta);
}

std::optional<double> HamiltonianReferenceEnergiesFo::get_es_exciation_enery([[maybe_unused]] unsigned n_k) const {
    return std::nullopt;
}

} // end of namespace model_monostar
