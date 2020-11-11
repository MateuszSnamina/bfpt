#include <so_hamiltonians/hamiltonian_reference_energies_af_fo.hpp>

#include <monostar_hamiltonians/hamiltonian_reference_energies_af_fm.hpp>

#include <cassert>

// #######################################################################
// ## helper: get_spin_reference_energies                               ##
// #######################################################################

namespace {
using namespace so_hamiltonians;

monostar_hamiltonians::HamiltonianReferenceEnergiesAf get_spin_reference_energies(
    HamiltonianParamsAfFo params,
    unsigned n_sites,
    double orbital_theta) {
    const monostar_hamiltonians::HamiltonianParamsAfFm spin_hamiltonian_params = params.average_out_orbitals_1(orbital_theta);
    monostar_hamiltonians::HamiltonianReferenceEnergiesAf spin_hamiltonian_reference_enrgies{n_sites, spin_hamiltonian_params};
    return spin_hamiltonian_reference_enrgies;
}

}  // end of namespace

// #######################################################################
// ## HamiltonianReferenceEnergiesFo                                    ##
// #######################################################################

namespace so_hamiltonians {

HamiltonianReferenceEnergiesAfFo::HamiltonianReferenceEnergiesAfFo(unsigned n_sites, const HamiltonianParamsAfFo& params, double orbital_theta)
    : HamiltonianReferenceEnergies(n_sites),
      _params(params),
      _orbital_theta(orbital_theta) {
}

std::optional<double> HamiltonianReferenceEnergiesAfFo::get_gs_energy() const {
    const auto spin_reference_energies = get_spin_reference_energies(_params, _n_sites, _orbital_theta);
    return spin_reference_energies.get_gs_energy();
}

std::optional<double> HamiltonianReferenceEnergiesAfFo::get_es_exciation_enery(unsigned n_k) const {
    const auto spin_reference_energies = get_spin_reference_energies(_params, _n_sites, _orbital_theta);
    return spin_reference_energies.get_es_exciation_enery(n_k);
}

}  // end of namespace so_hamiltonians
