#include <monostar_hamiltonians/hamiltonian_params_jkl01.hpp>

// #######################################################################
// ## HamiltonianParamsJkl01                                             ##
// #######################################################################

namespace monostar_hamiltonians {

HamiltonianParamsJkl01::Builder HamiltonianParamsJkl01::Builder::set_J(double J) {
    _J = J;
    return *this;
}

HamiltonianParamsJkl01::Builder HamiltonianParamsJkl01::Builder::set_K(double K) {
    _K = K;
    return *this;
}

HamiltonianParamsJkl01::Builder HamiltonianParamsJkl01::Builder::set_L(double L) {
    _L = L;
    return *this;
}

HamiltonianParamsJkl01::Builder HamiltonianParamsJkl01::Builder::set_J_0(double J_0) {
    _J_0 = J_0;
    return *this;
}

HamiltonianParamsJkl01::Builder HamiltonianParamsJkl01::Builder::set_K_0(double K_0) {
    _K_0 = K_0;
    return *this;
}

HamiltonianParamsJkl01::Builder HamiltonianParamsJkl01::Builder::set_L_0(double L_0) {
    _L_0 = L_0;
    return *this;
}

HamiltonianParamsJkl01::Builder HamiltonianParamsJkl01::Builder::set_J_1(double J_1) {
    _J_1 = J_1;
    return *this;
}

HamiltonianParamsJkl01::Builder HamiltonianParamsJkl01::Builder::set_K_1(double K_1) {
    _K_1 = K_1;
    return *this;
}

HamiltonianParamsJkl01::Builder HamiltonianParamsJkl01::Builder::set_L_1(double L_1) {
    _L_1 = L_1;
    return *this;
}

HamiltonianParamsJkl01 HamiltonianParamsJkl01::Builder::build() const {
    return HamiltonianParamsJkl01(_J, _K, _L,
                                  _J_0, _K_0, _L_0,
                                  _J_1, _K_1, _L_1);
}

HamiltonianParamsJkl01::HamiltonianParamsJkl01(
        double J, double K, double L,
        double J_0, double K_0, double L_0,
        double J_1, double K_1, double L_1)
    : _J(J), _K(K), _L(L),
      _J_0(J_0), _K_0(K_0), _L_0(L_0),
      _J_1(J_1), _K_1(K_1), _L_1(L_1) {
}

double HamiltonianParamsJkl01::get_J() const {
    return _J;
}

double HamiltonianParamsJkl01::get_K() const {
    return _K;
}

double HamiltonianParamsJkl01::get_L() const {
    return _L;
}

double HamiltonianParamsJkl01::get_J_0() const {
    return _J_0;
}

double HamiltonianParamsJkl01::get_K_0() const {
    return _K_0;
}

double HamiltonianParamsJkl01::get_L_0() const {
    return _L_0;
}

double HamiltonianParamsJkl01::get_J_1() const {
    return _J_1;
}

double HamiltonianParamsJkl01::get_K_1() const {
    return _K_1;
}

double HamiltonianParamsJkl01::get_L_1() const {
    return _L_1;
}

}  //end of namespace monostar_hamiltonians
