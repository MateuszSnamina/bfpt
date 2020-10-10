#pragma once

// #######################################################################
// ## HamiltonianParamsAfFm                                             ##
// #######################################################################

namespace monostar_hamiltonians {

class HamiltonianParamsAfFm {
public:
    double get_J_classical() const;
    double get_J_quantum() const;
    double get_B() const;
    class Builder {
    public:
        Builder set_J_classical(double);
        Builder set_J_quantum(double);
        Builder set_B(double);
        HamiltonianParamsAfFm build() const;
    private:
        double _J_classicall = 1.0;
        double _J_quantum = 1.0;
        double _B = 0.0;
    };
    friend HamiltonianParamsAfFm Builder::build() const;
    HamiltonianParamsAfFm() = default;
private:
    HamiltonianParamsAfFm(double J_classicall, double J_quantum, double B);
    double _J_classicall = 1.0;
    double _J_quantum = 1.0;
    double _B = 0.0;
};

} // end of namespace monostar_hamiltonians
