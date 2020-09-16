#ifndef MODEL_MONOSTAR_REFERENCE_ENERGIES_HPP
#define MODEL_MONOSTAR_REFERENCE_ENERGIES_HPP

#include<armadillo>

#include<cassert>

namespace model_monostar {

// #######################################################################
// ## ReferenceEnergies                                                 ##
// #######################################################################

class ReferenceEnergies {
public:
    ReferenceEnergies(unsigned n_sites) :
        _n_sites(n_sites) {
        assert(n_sites > 0);
    }
    virtual double get_gs_energy() const = 0;
    virtual double get_es_energy(unsigned n_k) const = 0;
    virtual ~ReferenceEnergies() = default;
protected:
    const unsigned _n_sites;
};

// #######################################################################
// ## ReferenceEnergiesFm                                               ##
// #######################################################################

class ReferenceEnergiesFm final : public ReferenceEnergies {
public:
    ReferenceEnergiesFm(unsigned n_sites, double J_classical, double J_quantum) :
        ReferenceEnergies(n_sites),
        _J_classical(J_classical),
        _J_quantum(J_quantum) {
        assert(_J_classical > 0);
    }
    double get_gs_energy() const override {
        return _J_classical * _n_sites * (- 0.25 );
    }
    double get_es_energy(unsigned n_k) const override {
        return _J_classical - _J_quantum * std::cos(2 * arma::datum::pi * n_k / _n_sites);
    }
private:
    const double _J_classical;
    const double _J_quantum;
};

// #######################################################################
// ## ReferenceEnergiesAf                                               ##
// #######################################################################

class ReferenceEnergiesAf final : public ReferenceEnergies {
public:
    ReferenceEnergiesAf(unsigned n_sites, double J) :
        ReferenceEnergies(n_sites),
        _J(J) {
    }
    double get_gs_energy() const override {
        return - _J * _n_sites * (std::log(2) - 0.25);
    }
    double get_es_energy(unsigned n_k) const override {
        return _J * arma::datum::pi/2 * std::abs(std::sin(2 * arma::datum::pi * n_k / _n_sites));
    }
private:
    const double _J;
};

}

#endif
