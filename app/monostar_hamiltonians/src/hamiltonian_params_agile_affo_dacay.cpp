#include <monostar_hamiltonians/hamiltonian_params_agile_affo_dacay.hpp>
#include <monostar_hamiltonians/hamiltonian_params_fo_site_matrices.hpp>

#include <armadillo>
#include <cmath>

//namespace  {

//arma::vec2 get_1e_psi_0() {
//    return arma::vec2{1.0, 0.0};
//}

//arma::vec4 get_2e_delta_psi(double phi) {
//    return arma::vec4{
//        0.0,
//        std::cos(phi) / std::sqrt(2),
//        std::cos(phi) / std::sqrt(2),
//        std::sin(phi)
//    };
//}

//arma::vec::fixed<16> get_4e_psi_0() {
//    const arma::vec2 _1e_psi_0 = get_1e_psi_0();
//    const arma::vec4 _2e_psi_0 = arma::kron(_1e_psi_0, _1e_psi_0);
//    const arma::vec8 _3e_psi_0 = arma::kron(_2e_psi_0, _1e_psi_0);
//    const arma::vec::fixed<16> _4e_psi_0 = arma::kron(_3e_psi_0, _1e_psi_0);
//    return _4e_psi_0;
//}

//arma::vec::fixed<16> get_4e_center_delta_psi(double phi) {
//    const arma::vec4 _2e_delta_psi = get_2e_delta_psi(phi);
//    const arma::vec2 _1e_psi_0 = get_1e_psi_0();
//    const arma::vec8 _3e_result = arma::kron(_1e_psi_0, _2e_delta_psi);
//    const arma::vec::fixed<16> _4e_result = arma::kron(_3e_result, _1e_psi_0);
//    return _4e_result;
//}

//arma::vec::fixed<16> get_4e_right_delta_psi(double phi) {
//    const arma::vec4 _2e_delta_psi = get_2e_delta_psi(phi);
//    const arma::vec2 _1e_psi_0 = get_1e_psi_0();
//    const arma::vec4 _2e_psi_0 = arma::kron(_1e_psi_0, _1e_psi_0);
//    const arma::vec::fixed<16> _4e_result = arma::kron(_2e_psi_0, _2e_delta_psi);
//    return _4e_result;
//}

//template<arma::uword N>
//double mean(arma::vec::fixed<N> state_1, arma::mat::fixed<N, N> op, arma::vec::fixed<N> state_2){
//    return arma::dot(state_1, op * state_2);
//}

//arma::mat::fixed<16, 16> operator_expander_2e_to_4e(arma::mat44 op_2e){
//    const E = arma::eye(2, 2);
//    const arma::mat88 _3e_op = arma::kron(E, op_2e);
//    const arma::mat::fixed<16, 16> _4e_op = arma::kron(op_3e, E);
//    return _4e_op;
//}

//double get_J() {
//    return mean (get_4e_psi_0(), get_4e_psi_0());
//}

////TODO finish!!

//}  // end of namespace

// #######################################################################
// ## dacay_hamiltonian_params_agile_affo                               ##
// #######################################################################

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
