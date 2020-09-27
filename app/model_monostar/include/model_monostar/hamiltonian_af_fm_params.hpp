#pragma once

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

inline
HamiltonianAfFmParams::HamiltonianAfFmParams(double J_classicall, double J_quantum, double B) :
    _J_classicall(J_classicall),
    _J_quantum(J_quantum),
    _B(B) {
}

inline
double HamiltonianAfFmParams::get_J_classical() const{
    return _J_classicall;
}

inline
double HamiltonianAfFmParams::get_J_quantum() const{
    return _J_quantum;
}

inline
double HamiltonianAfFmParams::get_B() const{
    return _B;
}

inline
HamiltonianAfFmParams::Builder HamiltonianAfFmParams::Builder::set_J_classical(double J_classicall) {
    _J_classicall = J_classicall;
    return *this;
}

inline
HamiltonianAfFmParams::Builder HamiltonianAfFmParams::Builder::set_J_quantum(double J_quantum) {
    _J_quantum = J_quantum;
    return *this;
}

inline
HamiltonianAfFmParams::Builder HamiltonianAfFmParams::Builder::set_B(double B) {
    _B = B;
    return *this;
}

inline
HamiltonianAfFmParams HamiltonianAfFmParams::Builder::build() const {
    return HamiltonianAfFmParams(_J_classicall, _J_quantum, _B);
}
