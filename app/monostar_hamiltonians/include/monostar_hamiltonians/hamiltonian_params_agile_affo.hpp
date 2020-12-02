#pragma once

#include <monostar_hamiltonians/hamiltonian_params_affo.hpp>

// #######################################################################
// ## AgileParams                                                       ##
// #######################################################################

namespace monostar_hamiltonians {

class AgileParams {
   public:
    class Builder {
       public:
        Builder set_eps(double);
        Builder set_phi(double);
        AgileParams build() const;

       private:
        double _eps = 0.0;
        double _phi = 0.0;
    };
    // Construction destruction:
    friend AgileParams Builder::build() const;
    // Getters:
    double get_eps() const;
    double get_phi() const;

   private:
    AgileParams(double eps, double phi);

   private:
    double _eps = 0.0;
    double _phi = 0.0;
};

}  // end of namespace monostar_hamiltonians

// #######################################################################
// ## HamiltonianParamsAgileAffo                                        ##
// #######################################################################

namespace monostar_hamiltonians {

class HamiltonianParamsAgileAffo {
   public:
    class Builder {
       public:
        Builder set_so_hamiltonian(HamiltonianParamsAffo);
        Builder set_aglie_params(AgileParams);
        HamiltonianParamsAgileAffo build() const;

       private:
        HamiltonianParamsAffo _so_hamiltonian = HamiltonianParamsAffo::Builder().build();
        AgileParams _aglie_params = AgileParams::Builder().build();
    };
    // Construction destruction:
    friend HamiltonianParamsAgileAffo Builder::build() const;
    // Getters:
    const HamiltonianParamsAffo& get_so_hamiltonian() const;
    const AgileParams& get_aglie_params() const;

   private:
    HamiltonianParamsAgileAffo(HamiltonianParamsAffo so_hamiltonian, AgileParams aglie_params);
    HamiltonianParamsAffo _so_hamiltonian;
    AgileParams _aglie_params;
};

}  // end of namespace monostar_hamiltonians
