#pragma once

#include<optional>

#include<cassert>

// #######################################################################
// ## HamiltonianReferenceEnergies                                      ##
// #######################################################################

namespace model_monostar {

class HamiltonianReferenceEnergies {
public:
    HamiltonianReferenceEnergies(unsigned n_sites) :
        _n_sites(n_sites) {
        assert(n_sites > 0);
    }
    virtual std::optional<double> get_gs_energy() const = 0;
    virtual std::optional<double> get_es_exciation_enery(unsigned n_k) const = 0;
    virtual std::optional<double> get_es_absolute_enery(unsigned n_k) const {
        if (const auto& gs_energy = get_gs_energy()) {
            if (const auto& es_exciation_enery = get_es_exciation_enery(n_k)) {
                return *gs_energy + *es_exciation_enery;
            }
        }
        return std::nullopt;
    }
    virtual ~HamiltonianReferenceEnergies() = default;
protected:
    const unsigned _n_sites;
};

}
