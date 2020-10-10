#include<monostar_hamiltonians/hamiltonian_reference_energies_af_fm.hpp>

#include<armadillo>

#include<optional>

#include<cassert>

// #######################################################################
// ## HamiltonianReferenceEnergiesFm                                    ##
// #######################################################################

namespace monostar_hamiltonians {

HamiltonianReferenceEnergiesFm::HamiltonianReferenceEnergiesFm(unsigned n_sites, double J_classical, double J_quantum, double B) :
    HamiltonianReferenceEnergies(n_sites),
    _J_classical(J_classical),
    _J_quantum(J_quantum),
    _B(B) {
    assert(_J_classical > 0);
}

HamiltonianReferenceEnergiesFm::HamiltonianReferenceEnergiesFm(unsigned n_sites, const HamiltonianParamsAfFm& params) :
    HamiltonianReferenceEnergies(n_sites),
    _J_classical(params.get_J_classical()),
    _J_quantum(params.get_J_quantum()),
    _B(params.get_B()) {
}

std::optional<double> HamiltonianReferenceEnergiesFm::get_gs_energy() const {
    return _n_sites * ((-0.25) * _J_classical + (-0.5) * _B);
}

std::optional<double> HamiltonianReferenceEnergiesFm::get_es_exciation_enery(unsigned n_k) const {
    return _J_classical + _B - _J_quantum * std::cos(2 * arma::datum::pi * n_k / _n_sites);
}

} // end of namespace monostar_hamiltonians

// #######################################################################
// ## HamiltonianReferenceEnergiesAf                                    ##
// #######################################################################

namespace monostar_hamiltonians {

HamiltonianReferenceEnergiesAf::HamiltonianReferenceEnergiesAf(unsigned n_sites, double J) :
    HamiltonianReferenceEnergies(n_sites),
    _J(J),
    _is_applicable(true){
}

HamiltonianReferenceEnergiesAf::HamiltonianReferenceEnergiesAf(unsigned n_sites, const HamiltonianParamsAfFm& params) :
    HamiltonianReferenceEnergies(n_sites),
    _J(params.get_J_classical()),
    _is_applicable(params.get_J_classical() == params.get_J_quantum() && params.get_B() == 0){
}

std::optional<double> HamiltonianReferenceEnergiesAf::get_gs_energy() const {
    if (_is_applicable) {
        return - _J * _n_sites * (std::log(2) - 0.25);
    } else {
        return std::nullopt;
    }
}

std::optional<double> HamiltonianReferenceEnergiesAf::get_es_exciation_enery(unsigned n_k) const {
    if (_is_applicable) {
        return _J * arma::datum::pi/2 * std::abs(std::sin(2 * arma::datum::pi * n_k / _n_sites));
    } else {
        return std::nullopt;
    }
}

} // end of namespace monostar_hamiltonians
