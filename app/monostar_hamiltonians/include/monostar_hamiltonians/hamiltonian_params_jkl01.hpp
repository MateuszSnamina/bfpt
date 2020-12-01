#pragma once

// #######################################################################
// ## HamiltonianParamsJkl01                                            ##
// #######################################################################

namespace monostar_hamiltonians {

class HamiltonianParamsJkl01 {
   public:
    class Builder {
       public:
        Builder set_J(double);
        Builder set_J_0(double);
        Builder set_J_1(double);
        Builder set_K(double);
        Builder set_K_0(double);
        Builder set_K_1(double);
        Builder set_L(double);
        Builder set_L_1(double);

        HamiltonianParamsJkl01 build() const;

       private:
        double _J = 1.0;
        double _J_0 = 0.0;
        double _J_1 = 0.0;
        double _K = 1.0;
        double _K_0 = 0.0;
        double _K_1 = 0.0;
        double _L = 0.0;
        double _L_1 = 0.0;
    };
    // Construction destruction:
    friend HamiltonianParamsJkl01 Builder::build() const;
    // Getters:
    double get_J() const;
    double get_K() const;
    double get_L() const;
    double get_J_0() const;
    double get_K_0() const;
    double get_L_0() const;
    double get_J_1() const;
    double get_K_1() const;
    double get_L_1() const;

   private:
    HamiltonianParamsJkl01(
            double J, double J_0, double J_1,
            double K, double K_0, double K_1,
            double L, double L_1);
    double _J;
    double _J_0;
    double _J_1;
    double _K;
    double _K_0;
    double _K_1;
    double _L;
    double _L_1;
};

}  // end of namespace monostar_hamiltonians
