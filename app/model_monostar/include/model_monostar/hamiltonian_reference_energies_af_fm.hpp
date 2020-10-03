#pragma once

#include<model_monostar/hamiltonian_reference_energies.hpp>
#include<model_monostar/hamiltonian_params_af_fm.hpp>

#include<armadillo>

#include<optional>

// #######################################################################
// ## HamiltonianReferenceEnergiesFm                                    ##
// #######################################################################

namespace model_monostar {

class HamiltonianReferenceEnergiesFm final : public HamiltonianReferenceEnergies {
public:
    HamiltonianReferenceEnergiesFm(unsigned n_sites, double J_classical, double J_quantum, double B);
    HamiltonianReferenceEnergiesFm(unsigned n_sites, const HamiltonianParamsAfFm& params);
    std::optional<double> get_gs_energy() const override;
    std::optional<double> get_es_exciation_enery(unsigned n_k) const override;
private:
    const double _J_classical;
    const double _J_quantum;
    const double _B;
};

} // end of namespace model_monostar

// #######################################################################
// ## HamiltonianReferenceEnergiesAf                                    ##
// #######################################################################

namespace model_monostar {

class HamiltonianReferenceEnergiesAf final : public HamiltonianReferenceEnergies {
public:
    HamiltonianReferenceEnergiesAf(unsigned n_sites, double J);
    HamiltonianReferenceEnergiesAf(unsigned n_sites, const HamiltonianParamsAfFm& params);
    std::optional<double> get_gs_energy() const override;
    std::optional<double> get_es_exciation_enery(unsigned n_k) const override;
private:
    const double _J;
    const bool _is_applicable;
};

} // end of namespace model_monostar
