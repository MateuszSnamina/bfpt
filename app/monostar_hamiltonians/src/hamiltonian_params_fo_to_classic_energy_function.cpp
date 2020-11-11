#include<monostar_hamiltonians/hamiltonian_params_fo_to_classic_energy_function.hpp>

// #######################################################################
// ## hamiltonian_params_fo_to_classic_energy_function                  ##
// #######################################################################

namespace monostar_hamiltonians {

AcosPlusBsinPlusCsqcosPlusZ hamiltonian_params_fo_to_classic_energy_function(HamiltonianParamsFo params) {
    // E(θ) = + A*cos(θ) + B*sin(θ)  + C*cos⁴(θ/2) + 2*D*cos²(θ/2)sin²(θ/2) + E*sin⁴(θ/2) + F
    //      = + A*cos(θ) + B*sin(θ)  + C[½+½cos(θ)]² + 2*D[½*sin(θ)]² + E[½-½cos(θ)]² + F
    //      = + A*cos(θ) + B*sin(θ)  + ¼C[1+cos(θ)]² + ½D[sin(θ)]² + ¼E[1-cos(θ)]² + F
    //      = + A*cos(θ) + B*sin(θ)  + (¼C+¼E) + ½(C-E)*cos(θ) + ¼(C+E)*cos²(θ) + ½D - ½D*cos²(θ) + F
    //      =  [A+½(C-E)]*cos(θ) - B*sin(θ) + ¼(C+E-2D)*cos²(θ) + ¼(C+E+2D) + F
    // where: A ≡ +tau_z_coef,
    //        B ≡ -tau_minus_coef,
    //        C ≡ +Pzz_coef,
    //        D ≡ +Pxz_coef,
    //        E ≡ +Pxx_coef.
    //        E ≡ +free_coef.
    const double cos_coef = params.get_tau_z_coef() + 0.5 * (+params.get_Pzz_coef() - params.get_Pxx_coef());
    const double sin_coef = -params.get_tau_minus_coef();
    const double sqcos_coef = 0.25 * (params.get_Pzz_coef() + params.get_Pxx_coef() - 2 * params.get_Pxz_coef());
    const double free_coef = 0.25 * (params.get_Pzz_coef() + params.get_Pxx_coef() + 2 * params.get_Pxz_coef()) + params.get_free_coef();
    return AcosPlusBsinPlusCsqcosPlusZ::Builder()
        .set_cos_coef(cos_coef)
        .set_sin_coef(sin_coef)
        .set_sqcos_coef(sqcos_coef)
        .set_free_coef(free_coef)
        .build();
}

}  // end of namespace monostar_hamiltonians
