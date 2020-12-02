#include <monostar_hamiltonians/hamiltonian_params_agile_affo.hpp>

// #######################################################################
// ## AgileParams                                                       ##
// #######################################################################

namespace monostar_hamiltonians {

AgileParams::Builder AgileParams::Builder::set_eps(double eps) {
    _eps = eps;
    return *this;
}

AgileParams::Builder AgileParams::Builder::set_phi(double phi) {
    _phi = phi;
    return *this;
}

AgileParams AgileParams::Builder::build() const {
    return AgileParams(_eps, _phi);
}

AgileParams::AgileParams(double eps, double phi)
    : _eps(eps),
      _phi(phi) {
}

double AgileParams::get_eps() const {
    return _eps;
}

double AgileParams::get_phi() const {
    return _phi;
}

}  //end of namespace monostar_hamiltonians

// #######################################################################
// ## HamiltonianParamsAgileAffo                                        ##
// #######################################################################

namespace monostar_hamiltonians {

HamiltonianParamsAgileAffo::Builder HamiltonianParamsAgileAffo::Builder::set_so_hamiltonian(HamiltonianParamsAffo so_hamiltonian) {
    _so_hamiltonian = so_hamiltonian;
    return *this;
}

HamiltonianParamsAgileAffo::Builder HamiltonianParamsAgileAffo::Builder::set_aglie_params(AgileParams aglie_params) {
    _aglie_params = aglie_params;
    return *this;
}

HamiltonianParamsAgileAffo HamiltonianParamsAgileAffo::Builder::build() const {
    return HamiltonianParamsAgileAffo(_so_hamiltonian, _aglie_params);
}

HamiltonianParamsAgileAffo::HamiltonianParamsAgileAffo(HamiltonianParamsAffo so_hamiltonian, AgileParams aglie_params)
    : _so_hamiltonian(so_hamiltonian),
      _aglie_params(aglie_params) {
}

const HamiltonianParamsAffo& HamiltonianParamsAgileAffo::get_so_hamiltonian() const {
    return _so_hamiltonian;
}

const AgileParams& HamiltonianParamsAgileAffo::get_aglie_params() const {
    return _aglie_params;
}

}  //end of namespace monostar_hamiltonians
