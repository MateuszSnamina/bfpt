#include <monostar_hamiltonians/hamiltonian_params_af_fm.hpp>

// #######################################################################
// ## HamiltonianParamsAfFm                                             ##
// #######################################################################

namespace monostar_hamiltonians {

HamiltonianParamsAfFm::Builder HamiltonianParamsAfFm::Builder::set_J_classical(double J_classicall) {
    _J_classicall = J_classicall;
    return *this;
}

HamiltonianParamsAfFm::Builder HamiltonianParamsAfFm::Builder::set_J_quantum(double J_quantum) {
    _J_quantum = J_quantum;
    return *this;
}

HamiltonianParamsAfFm::Builder HamiltonianParamsAfFm::Builder::set_B(double B) {
    _B = B;
    return *this;
}

HamiltonianParamsAfFm::Builder HamiltonianParamsAfFm::Builder::set_free(double free) {
    _free = free;
    return *this;
}

HamiltonianParamsAfFm HamiltonianParamsAfFm::Builder::build() const {
    return HamiltonianParamsAfFm(_J_classicall, _J_quantum, _B, _free);
}

HamiltonianParamsAfFm::HamiltonianParamsAfFm(double J_classicall, double J_quantum, double B, double free)
    : _J_classicall(J_classicall),
      _J_quantum(J_quantum),
      _B(B),
      _free(free) {
}

double HamiltonianParamsAfFm::get_J_classical() const {
    return _J_classicall;
}

double HamiltonianParamsAfFm::get_J_quantum() const {
    return _J_quantum;
}

double HamiltonianParamsAfFm::get_B() const {
    return _B;
}

double HamiltonianParamsAfFm::get_free() const {
    return _free;
}

}  //end of namespace monostar_hamiltonians
