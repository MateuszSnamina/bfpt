#pragma once

// #######################################################################
// ## HamiltonianAfFmParams                                             ##
// #######################################################################

class HamiltonianAfFmParams {
public:
    double get_J_classical() const;
    double get_J_quantum() const;
    double get_B() const;
    class Builder {
    public:
        Builder set_J_classical(double);
        Builder set_J_quantum(double);
        Builder set_B(double);
        HamiltonianAfFmParams build() const;
    private:
        double _J_classicall = 1.0;
        double _J_quantum = 1.0;
        double _B = 0.0;
    };
    friend HamiltonianAfFmParams Builder::build() const;
    HamiltonianAfFmParams() = default;
private:
    HamiltonianAfFmParams(double J_classicall, double J_quantum, double B);
    double _J_classicall = 1.0;
    double _J_quantum = 1.0;
    double _B = 0.0;
};
