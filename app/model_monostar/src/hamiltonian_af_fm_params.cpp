#include<model_monostar/hamiltonian_af_fm_params.hpp>

HamiltonianAfFmParams::Builder HamiltonianAfFmParams::Builder::set_J_classical(double J_classicall) {
    _J_classicall = J_classicall;
    return *this;
}

HamiltonianAfFmParams::Builder HamiltonianAfFmParams::Builder::set_J_quantum(double J_quantum) {
    _J_quantum = J_quantum;
    return *this;
}

HamiltonianAfFmParams::Builder HamiltonianAfFmParams::Builder::set_B(double B) {
    _B = B;
    return *this;
}

HamiltonianAfFmParams HamiltonianAfFmParams::Builder::build() const {
    return HamiltonianAfFmParams(_J_classicall, _J_quantum, _B);
}

HamiltonianAfFmParams::HamiltonianAfFmParams(double J_classicall, double J_quantum, double B) :
    _J_classicall(J_classicall),
    _J_quantum(J_quantum),
    _B(B) {
}

double HamiltonianAfFmParams::get_J_classical() const{
    return _J_classicall;
}

double HamiltonianAfFmParams::get_J_quantum() const{
    return _J_quantum;
}

double HamiltonianAfFmParams::get_B() const{
    return _B;
}
