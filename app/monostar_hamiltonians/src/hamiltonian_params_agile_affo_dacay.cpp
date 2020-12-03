#include<monostar_hamiltonians/hamiltonian_params_agile_affo_dacay.hpp>
#include<monostar_hamiltonians/hamiltonian_params_fo_site_matrices.hpp>

namespace monostar_hamiltonians {

HamiltonianParamsJkl01 dacay_hamiltonian_params_agile_affo(
        HamiltonianParamsAgileAffo hamiltonian_params_agile_affo,
        double orbital_theta) {

    const arma::cx_mat44 P_zz = TwoSitesOrbitalMatrices::get_P_zz_in_ge_basis(orbital_theta);
    const arma::cx_mat44 P_zx_sum_P_xz = TwoSitesOrbitalMatrices::get_P_zx_sum_P_xz_in_ge_basis(orbital_theta);
    const arma::cx_mat44 P_xx = TwoSitesOrbitalMatrices::get_P_xx_in_ge_basis(orbital_theta);

    const auto so_hamiltonian_pure = hamiltonian_params_agile_affo.get_so_hamiltonian().get_hamiltonian_params_fo();
    const arma::cx_mat44 orbital_operator_pure =
            so_hamiltonian_pure.get_Pzz_coef() * P_zz +
            so_hamiltonian_pure.get_Pxz_coef() * P_zx_sum_P_xz +
            so_hamiltonian_pure.get_Pxx_coef() * P_xx;
    const auto so_hamiltonian_prefactored_ss = hamiltonian_params_agile_affo.get_so_hamiltonian().get_hamiltonian_params_ss_fo();
    const arma::cx_mat44 orbital_operator_prefactored_ss =
            so_hamiltonian_prefactored_ss.get_Pzz_coef() * P_zz +
            so_hamiltonian_prefactored_ss.get_Pxz_coef() * P_zx_sum_P_xz +
            so_hamiltonian_prefactored_ss.get_Pxx_coef() * P_xx;


    const double J = hamiltonian_params_agile_affo.get_so_hamiltonian().get_ss_coef();
    const double J_0 = 0.0;
    const double J_1 = 0.0;
    const double K = 0.0;
    const double K_0 = 0.0;
    const double K_1 = 0.0;
    const double L = 0.0;
    const double L_1 = 0.0;

    // TODO: FINISH !!!

    return monostar_hamiltonians::HamiltonianParamsJkl01::Builder()
            .set_J(J)
            .set_J_0(J_0)
            .set_J_1(J_1)
            .set_K(K)
            .set_K_0(K_0)
            .set_K_1(K_1)
            .set_L(L)
            .set_L_1(L_1)
            .build();
}

}  // end of namespace monostar_hamiltonians
