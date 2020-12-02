#include <monostar_hamiltonians/hamiltonian_reference_energies_jkl01.hpp>

#include <monostar_hamiltonians/hamiltonian_params_af_fm.hpp>
#include <monostar_hamiltonians/hamiltonian_reference_energies_af_fm.hpp>

#include <cassert>

// #######################################################################
// ## helpers                                                           ##
// #######################################################################

namespace {

using namespace monostar_hamiltonians;
HamiltonianParamsAfFm hamiltonian_params_jkl01_to_Hamiltonian_params_af(HamiltonianParamsJkl01 params_jkl01) {
    return HamiltonianParamsAfFm::Builder()
        .set_J_classical(params_jkl01.get_J())
        .set_J_quantum(params_jkl01.get_J())
        .set_free(params_jkl01.get_K() + params_jkl01.get_L())
        .build();
}

}  // end of namespace

// #######################################################################
// ## HamiltonianReferenceEnergiesJkl01                                 ##
// #######################################################################

namespace monostar_hamiltonians {

HamiltonianReferenceEnergiesJkl01::HamiltonianReferenceEnergiesJkl01(unsigned n_sites, const HamiltonianParamsJkl01& params)
    : HamiltonianReferenceEnergies(n_sites),
      _params(params) {
}

std::optional<double> HamiltonianReferenceEnergiesJkl01::get_gs_energy() const {
    const HamiltonianParamsAfFm hamiltonian_params_af = hamiltonian_params_jkl01_to_Hamiltonian_params_af(_params);
    const HamiltonianReferenceEnergiesAf hamiltonian_reference_energies(_n_sites, hamiltonian_params_af);
    return hamiltonian_reference_energies.get_gs_energy();
}

std::optional<double> HamiltonianReferenceEnergiesJkl01::get_es_exciation_enery(unsigned n_k) const {
    const HamiltonianParamsAfFm hamiltonian_params_af = hamiltonian_params_jkl01_to_Hamiltonian_params_af(_params);
    const HamiltonianReferenceEnergiesAf hamiltonian_reference_energies(_n_sites, hamiltonian_params_af);
    return hamiltonian_reference_energies.get_es_exciation_enery(n_k);
}

}  // end of namespace monostar_hamiltonians
