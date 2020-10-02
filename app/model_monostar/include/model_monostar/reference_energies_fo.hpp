#pragma once

#include<model_monostar/reference_energies.hpp>
#include<model_monostar/hamiltonian_fo_params.hpp>

#include<cassert>

//// #######################################################################
//// ## ReferenceEnergiesFm                                               ##
//// #######################################################################

//namespace model_monostar {

//class ReferenceEnergiesFm final : public ReferenceEnergies {
//public:
//    ReferenceEnergiesFm(unsigned n_sites, double J_classical, double J_quantum, double B) :
//        ReferenceEnergies(n_sites),
//        _J_classical(J_classical),
//        _J_quantum(J_quantum),
//        _B(B) {
//        assert(_J_classical > 0);
//    }
//    std::optional<double> get_gs_energy() const override {
//        return _n_sites * ((-0.25) * _J_classical + (-0.5) * _B);
//    }
//    std::optional<double> get_es_exciation_enery(unsigned n_k) const override {
//        return _J_classical + _B - _J_quantum * std::cos(2 * arma::datum::pi * n_k / _n_sites);
//    }
//private:
//    const double _J_classical;
//    const double _J_quantum;
//    const double _B;
//};

//} // end of namespace model_monostar

//// #######################################################################
//// ## ReferenceEnergiesAf                                               ##
//// #######################################################################

//namespace model_monostar {

//class ReferenceEnergiesAf final : public ReferenceEnergies {
//public:
//    ReferenceEnergiesAf(unsigned n_sites, double J) :
//        ReferenceEnergies(n_sites),
//        _J(J) {
//    }
//    std::optional<double> get_gs_energy() const override {
//        return - _J * _n_sites * (std::log(2) - 0.25);
//    }
//    std::optional<double> get_es_exciation_enery(unsigned n_k) const override {
//        return _J * arma::datum::pi/2 * std::abs(std::sin(2 * arma::datum::pi * n_k / _n_sites));
//    }
//private:
//    const double _J;
//};

//} // end of namespace model_monostar
