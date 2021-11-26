#include <monostar_hamiltonians/hamiltonian_reference_energies_af_fm.hpp>

#include <armadillo>

#include <optional>

#include <cassert>

// #######################################################################
// ## HamiltonianReferenceEnergiesFm                                    ##
// #######################################################################

namespace monostar_hamiltonians {

HamiltonianReferenceEnergiesFm::HamiltonianReferenceEnergiesFm(unsigned n_sites,
                                                               double J_classical, double J_quantum,
                                                               double J_nnn_classical, double J_nnn_quantum,
                                                               double B, double free)
    : HamiltonianReferenceEnergies(n_sites),
      _J_classical(J_classical),
      _J_quantum(J_quantum),
      _J_nnn_classical(J_nnn_classical),
      _J_nnn_quantum(J_nnn_quantum),
      _B(B),
      _free(free) {
    assert(_J_classical > 0);
}

HamiltonianReferenceEnergiesFm::HamiltonianReferenceEnergiesFm(unsigned n_sites, const HamiltonianParamsAfFm& params)
    : HamiltonianReferenceEnergies(n_sites),
      _J_classical(params.get_J_classical()),
      _J_quantum(params.get_J_quantum()),
      _J_nnn_classical(params.get_J_nnn_classical()),
      _J_nnn_quantum(params.get_J_nnn_quantum()),
      _B(params.get_B()),
      _free(params.get_free()) {
}

std::optional<double> HamiltonianReferenceEnergiesFm::get_gs_energy() const {
    // TODO check: |_J_classical| > |_J_nnn_classical| condition!
    return _n_sites * (+(-0.25) * _J_classical + (+0.25) * _J_nnn_classical + (-0.5) * _B + _free);
}

std::optional<double> HamiltonianReferenceEnergiesFm::get_es_exciation_enery(unsigned n_k) const {
    return +_J_classical - _J_nnn_classical + _B - _J_quantum * std::cos(2 * arma::datum::pi * n_k / _n_sites) + _J_nnn_quantum * std::cos(2 * 2 * arma::datum::pi * n_k / _n_sites);
}

}  // end of namespace monostar_hamiltonians

// #######################################################################
// ## HamiltonianReferenceEnergiesAf                                    ##
// #######################################################################

namespace monostar_hamiltonians {

HamiltonianReferenceEnergiesAf::HamiltonianReferenceEnergiesAf(unsigned n_sites, double J, double free)
    : HamiltonianReferenceEnergies(n_sites),
      _J(J),
      _free(free),
      _is_applicable(true) {
}

HamiltonianReferenceEnergiesAf::HamiltonianReferenceEnergiesAf(unsigned n_sites, const HamiltonianParamsAfFm& params)
    : HamiltonianReferenceEnergies(n_sites),
      _J(params.get_J_classical()),
      _free(params.get_free()),
      _is_applicable(params.get_J_classical() == params.get_J_quantum() &&
                     params.get_B() == 0 &&
                     params.get_J_nnn_classical() == 0 &&
                     params.get_J_nnn_quantum() == 0) {
}

std::optional<double> HamiltonianReferenceEnergiesAf::get_gs_energy() const {
    if (_is_applicable) {
        return -_J * _n_sites * (std::log(2) - 0.25) + _n_sites * _free;
    } else {
        return std::nullopt;
    }
}

std::optional<double> HamiltonianReferenceEnergiesAf::get_es_exciation_enery(unsigned n_k) const {
    if (_is_applicable) {
        return _J * arma::datum::pi / 2 * std::abs(std::sin(2 * arma::datum::pi * n_k / _n_sites));
    } else {
        return std::nullopt;
    }
}

}  // end of namespace monostar_hamiltonians
