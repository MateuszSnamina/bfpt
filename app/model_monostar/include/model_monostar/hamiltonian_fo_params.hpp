#pragma once

class HamiltonianFoParams {
public:
    double get_Pdelta_coef() const;
    double get_Pxx_coef() const;
    double get_Pxz_coef() const;
    double get_Pzz_coef() const;
    class Builder {
    public:
        Builder set_Pdelta_coef(double);
        Builder set_Pxx_coef(double);
        Builder set_Pxz_coef(double);
        Builder set_Pzz_coef(double);
        HamiltonianFoParams build() const;
    private:
        double _Pdelta_coef = 1.0;
        double _Pxx_coef = 1.0;
        double _Pxz_coef = 0.0;
        double _Pzz_coef = 0.0;
    };
    friend HamiltonianFoParams Builder::build() const;
    HamiltonianFoParams() = default;
private:
    HamiltonianFoParams(double Bdelta, double Pxx_coef, double Pxz_coef, double Pzz_coef);
    double _Pdelta_coef = 1.0;
    double _Pxx_coef = 1.0;
    double _Pxz_coef = 0.0;
    double _Pzz_coef = 0.0;

};

inline
HamiltonianFoParams::HamiltonianFoParams(double Pdelta_coef, double Pxx_coef, double Pxz_coef, double Pzz_coef) :
    _Pdelta_coef(Pdelta_coef),
    _Pxx_coef(Pxx_coef),
    _Pxz_coef(Pxz_coef),
    _Pzz_coef(Pzz_coef){
}

inline
double HamiltonianFoParams::get_Pdelta_coef() const{
    return _Pdelta_coef;
}

inline
double HamiltonianFoParams::get_Pxx_coef() const{
    return _Pxx_coef;
}

inline
double HamiltonianFoParams::get_Pxz_coef() const{
    return _Pxz_coef;
}

inline
double HamiltonianFoParams::get_Pzz_coef() const{
    return _Pzz_coef;
}


inline
HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_Pdelta_coef(double Pdelta_coef) {
    _Pdelta_coef = Pdelta_coef;
    return *this;
}

inline
HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_Pxx_coef(double Pxx_coef) {
    _Pxx_coef = Pxx_coef;
    return *this;
}

inline
HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_Pxz_coef(double Pxz_coef) {
    _Pxz_coef = Pxz_coef;
    return *this;
}

inline
HamiltonianFoParams::Builder HamiltonianFoParams::Builder::set_Pzz_coef(double Pzz_coef) {
    _Pzz_coef = Pzz_coef;
    return *this;
}

inline
HamiltonianFoParams HamiltonianFoParams::Builder::build() const {
    return HamiltonianFoParams(_Pdelta_coef, _Pxx_coef, _Pxz_coef, _Pzz_coef);
}
