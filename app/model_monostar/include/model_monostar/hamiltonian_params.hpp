#pragma once

class HamiltonianParams {
public:
    double get_J_classical() const;
    double get_J_quantum() const;
    double get_B() const;
    class Builder {
    public:
        Builder set_J_classical(double);
        Builder set_J_quantum(double);
        Builder set_B(double);
        HamiltonianParams build() const;
    private:
        double _J_classicall = 1.0;
        double _J_quantum = 1.0;
        double _B = 0.0;
    };
    friend HamiltonianParams Builder::build() const;
    HamiltonianParams() = default;
private:
    HamiltonianParams(double J_classicall, double J_quantum, double B);
    double _J_classicall = 1.0;
    double _J_quantum = 1.0;
    double _B = 0.0;
};

inline
HamiltonianParams::HamiltonianParams(double J_classicall, double J_quantum, double B) :
    _J_classicall(J_classicall),
    _J_quantum(J_quantum),
    _B(B) {
}

inline
double HamiltonianParams::get_J_classical() const{
    return _J_classicall;
}

inline
double HamiltonianParams::get_J_quantum() const{
    return _J_quantum;
}

inline
double HamiltonianParams::get_B() const{
    return _B;
}

inline
HamiltonianParams::Builder HamiltonianParams::Builder::set_J_classical(double J_classicall) {
    _J_classicall = J_classicall;
    return *this;
}

inline
HamiltonianParams::Builder HamiltonianParams::Builder::set_J_quantum(double J_quantum) {
    _J_quantum = J_quantum;
    return *this;
}

inline
HamiltonianParams::Builder HamiltonianParams::Builder::set_B(double B) {
    _B = B;
    return *this;
}

inline
HamiltonianParams HamiltonianParams::Builder::build() const {
    return HamiltonianParams(_J_classicall, _J_quantum, _B);
}
