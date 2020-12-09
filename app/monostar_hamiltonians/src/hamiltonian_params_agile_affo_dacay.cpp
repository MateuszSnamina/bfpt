#include <monostar_hamiltonians/hamiltonian_params_agile_affo_dacay.hpp>
#include <monostar_hamiltonians/hamiltonian_params_fo_site_matrices.hpp>

#include <armadillo>
#include <cmath>

namespace {

/*
 * Implementation note:
 * All of the calculations below for orbital operator
 * matrix elements are performed in ge-basis (not in zx-basis)
 */

arma::cx_vec2 get_1e_gs() {
    return arma::cx_vec2{1.0, 0.0};
}

arma::cx_vec4 get_2e_psi(double phi) {
    return arma::cx_vec4{
        0.0,
        std::cos(phi) / std::sqrt(2),
        std::cos(phi) / std::sqrt(2),
        std::sin(phi)};
}

arma::cx_vec8 get_state_gs_x_gs_x_gs() {
    // Get brick:
    const arma::cx_vec2 _1e_gs = get_1e_gs();
    // Build tensor product:
    const arma::cx_vec4 gs_x_gs = arma::kron(_1e_gs, _1e_gs);
    const arma::cx_vec8 gs_x_gs_x_gs = arma::kron(gs_x_gs, _1e_gs);
    return gs_x_gs_x_gs;
}

arma::cx_vec8 get_state_gs_x_psi(double phi) {
    // Get bricks:
    const arma::cx_vec2 _1e_gs = get_1e_gs();
    const arma::cx_vec4 _2e_psi = get_2e_psi(phi);
    // Build tensor product:
    const arma::cx_vec8 gs_x_psi = arma::kron(_1e_gs, _2e_psi);
    return gs_x_psi;
}

arma::cx_vec::fixed<16> get_state_gs_x_gs_x_gs_x_gs() {
    // Get brick:
    const arma::cx_vec2 _1e_gs = get_1e_gs();
    // Build tensor product:
    const arma::cx_vec4 gs_x_gs = arma::kron(_1e_gs, _1e_gs);
    const arma::cx_vec8 gs_x_gs_x_gs = arma::kron(gs_x_gs, _1e_gs);
    const arma::cx_vec::fixed<16> gs_x_gs_x_gs_x_gs = arma::kron(gs_x_gs_x_gs, _1e_gs);
    return gs_x_gs_x_gs_x_gs;
}

arma::cx_vec::fixed<16> get_state_gs_x_psi_x_gs(double phi) {
    // Get bricks:
    const arma::cx_vec2 _1e_gs = get_1e_gs();
    const arma::cx_vec4 _2e_psi = get_2e_psi(phi);
    // Build tensor product:
    const arma::cx_vec8 gs_x_psi = arma::kron(_1e_gs, _2e_psi);
    const arma::cx_vec::fixed<16> gs_x_psi_x_gs = arma::kron(gs_x_psi, _1e_gs);
    return gs_x_psi_x_gs;
}

arma::cx_vec::fixed<16> get_state_gs_x_gs_x_psi(double phi) {
    // Get bricks:
    const arma::cx_vec2 _1e_gs = get_1e_gs();
    const arma::cx_vec4 _2e_psi = get_2e_psi(phi);
    // Build tensor product:
    const arma::cx_vec4 gs_x_gs = arma::kron(_1e_gs, _1e_gs);
    const arma::cx_vec::fixed<16> gs_x_gs_x_psi = arma::kron(gs_x_gs, _2e_psi);
    return gs_x_gs_x_psi;
}

template <arma::uword N>
double mean(arma::cx_vec::fixed<N> state_1, arma::cx_mat::fixed<N, N> op, arma::cx_vec::fixed<N> state_2) {
    std::complex<double> result = arma::dot(state_1, op * state_2);
    if (std::abs(std::imag(result)) > 1e-7) {
        throw std::runtime_error("Internal error");  //TODO make error handling nice.
    }
    return std::real(result);
}

template <arma::uword N>
arma::cx_mat::fixed<4u * N, 4u * N> get_id_x_op_x_id(arma::cx_mat::fixed<N, N> op) {
    const arma::cx_mat22 id_1e = arma::cx_mat22(arma::eye(2, 2), arma::zeros(2, 2));
    const arma::cx_mat::fixed<2u * N, 2u * N> id_x_op = arma::kron(id_1e, op);
    const arma::cx_mat::fixed<4u * N, 4u * N> id_x_op_x_id = arma::kron(id_x_op, id_1e);
    return id_x_op_x_id;
}

double get_J(const arma::cx_mat44& orbital_operator_2e_prefactored_ss) {
    return mean<16>(get_state_gs_x_gs_x_gs_x_gs(),
                    get_id_x_op_x_id<4>(orbital_operator_2e_prefactored_ss),
                    get_state_gs_x_gs_x_gs_x_gs());
}

double get_J_0(const arma::cx_mat44& orbital_operator_2e_prefactored_ss,
               const double eps, const double phi) {
    return +2 * eps * mean<16>(get_state_gs_x_gs_x_gs_x_gs(), get_id_x_op_x_id<4>(orbital_operator_2e_prefactored_ss), get_state_gs_x_psi_x_gs(phi)) + eps * eps * mean<16>(get_state_gs_x_psi_x_gs(phi), get_id_x_op_x_id<4>(orbital_operator_2e_prefactored_ss), get_state_gs_x_psi_x_gs(phi)) - eps * eps * mean<16>(get_state_gs_x_gs_x_gs_x_gs(), get_id_x_op_x_id<4>(orbital_operator_2e_prefactored_ss), get_state_gs_x_gs_x_gs_x_gs());
}

double get_J_1(const arma::cx_mat44& orbital_operator_2e_prefactored_ss,
               const double eps, const double phi) {
    return 2 * eps * mean<16>(get_state_gs_x_gs_x_gs_x_gs(), get_id_x_op_x_id<4>(orbital_operator_2e_prefactored_ss), get_state_gs_x_gs_x_psi(phi)) + eps * eps * mean<16>(get_state_gs_x_gs_x_psi(phi), get_id_x_op_x_id<4>(orbital_operator_2e_prefactored_ss), get_state_gs_x_gs_x_psi(phi)) - eps * eps * mean<16>(get_state_gs_x_gs_x_gs_x_gs(), get_id_x_op_x_id<4>(orbital_operator_2e_prefactored_ss), get_state_gs_x_gs_x_gs_x_gs());
}

double get_K(const arma::cx_mat44& orbital_operator_2e_pure) {
    return mean<16>(get_state_gs_x_gs_x_gs_x_gs(),
                    get_id_x_op_x_id<4>(orbital_operator_2e_pure),
                    get_state_gs_x_gs_x_gs_x_gs());
}

double get_K_0(const arma::cx_mat44& orbital_operator_2e_pure,
               const double eps, const double phi) {
    return 2 * eps * mean<16>(get_state_gs_x_gs_x_gs_x_gs(), get_id_x_op_x_id<4>(orbital_operator_2e_pure), get_state_gs_x_psi_x_gs(phi)) + eps * eps * mean<16>(get_state_gs_x_psi_x_gs(phi), get_id_x_op_x_id<4>(orbital_operator_2e_pure), get_state_gs_x_psi_x_gs(phi)) - eps * eps * mean<16>(get_state_gs_x_gs_x_gs_x_gs(), get_id_x_op_x_id<4>(orbital_operator_2e_pure), get_state_gs_x_gs_x_gs_x_gs());
}

double get_K_1(const arma::cx_mat44& orbital_operator_2e_pure,
               const double eps, const double phi) {
    return 2 * eps * mean<16>(get_state_gs_x_gs_x_gs_x_gs(), get_id_x_op_x_id<4>(orbital_operator_2e_pure), get_state_gs_x_gs_x_psi(phi)) + eps * eps * mean<16>(get_state_gs_x_gs_x_psi(phi), get_id_x_op_x_id<4>(orbital_operator_2e_pure), get_state_gs_x_gs_x_psi(phi)) - eps * eps * mean<16>(get_state_gs_x_gs_x_gs_x_gs(), get_id_x_op_x_id<4>(orbital_operator_2e_pure), get_state_gs_x_gs_x_gs_x_gs());
}

double get_L(const arma::cx_mat22& orbital_operator_1e_pure) {
    return mean<8>(get_state_gs_x_gs_x_gs(),
                   get_id_x_op_x_id<2>(orbital_operator_1e_pure),
                   get_state_gs_x_gs_x_gs());
}

double get_L_1(const arma::cx_mat22& orbital_operator_1e_pure,
               const double eps, const double phi) {
    return 2 * eps * mean<8>(get_state_gs_x_gs_x_gs(), get_id_x_op_x_id<2>(orbital_operator_1e_pure), get_state_gs_x_psi(phi)) + eps * eps * mean<8>(get_state_gs_x_psi(phi), get_id_x_op_x_id<2>(orbital_operator_1e_pure), get_state_gs_x_psi(phi)) - eps * eps * mean<8>(get_state_gs_x_gs_x_gs(), get_id_x_op_x_id<2>(orbital_operator_1e_pure), get_state_gs_x_gs_x_gs());
}

}  // end of namespace

// #######################################################################
// ## dacay_hamiltonian_params_agile_affo                               ##
// #######################################################################

namespace monostar_hamiltonians {

HamiltonianParamsJkl01 dacay_hamiltonian_params_agile_affo(
    HamiltonianParamsAgileAffo hamiltonian_params_agile_affo,
    double orbital_theta) {
    // Helper matrices:
    const arma::cx_mat22& tau_minus = OneSiteOrbitalMatrices::get_tau_minus_in_ge_basis(orbital_theta);
    const arma::cx_mat22& tau_z = OneSiteOrbitalMatrices::get_tau_z_in_ge_basis(orbital_theta);
    const arma::cx_mat44& P_zz = TwoSitesOrbitalMatrices::get_P_zz_in_ge_basis(orbital_theta);
    const arma::cx_mat44& P_zx_sum_P_xz = TwoSitesOrbitalMatrices::get_P_zx_sum_P_xz_in_ge_basis(orbital_theta);
    const arma::cx_mat44& P_xx = TwoSitesOrbitalMatrices::get_P_xx_in_ge_basis(orbital_theta);
    // Refs for hamiltonian_params_agile_affo items:
    const auto& so_hamiltonian_prefactored_ss = hamiltonian_params_agile_affo.get_so_hamiltonian().get_hamiltonian_params_ss_fo();
    const auto& so_hamiltonian_pure = hamiltonian_params_agile_affo.get_so_hamiltonian().get_hamiltonian_params_fo();
    const auto& eps = hamiltonian_params_agile_affo.get_aglie_params().get_eps();
    const auto& phi = hamiltonian_params_agile_affo.get_aglie_params().get_phi();
    // Helper orbital operators:
    const arma::cx_mat44 orbital_operator_2e_prefactored_ss =
        so_hamiltonian_prefactored_ss.get_Pzz_coef() * P_zz +
        so_hamiltonian_prefactored_ss.get_Pxz_coef() * P_zx_sum_P_xz +
        so_hamiltonian_prefactored_ss.get_Pxx_coef() * P_xx;
    const arma::cx_mat44 orbital_operator_2e_pure =
        so_hamiltonian_pure.get_Pzz_coef() * P_zz +
        so_hamiltonian_pure.get_Pxz_coef() * P_zx_sum_P_xz +
        so_hamiltonian_pure.get_Pxx_coef() * P_xx;
    const arma::cx_mat22 orbital_operator_1e_pure =
        so_hamiltonian_pure.get_tau_minus_coef() * tau_minus +
        so_hamiltonian_pure.get_tau_z_coef() * tau_z;
    // Jkl01 Hamiltonian params:
    const double J = get_J(orbital_operator_2e_prefactored_ss) + hamiltonian_params_agile_affo.get_so_hamiltonian().get_ss_coef();
    const double J_0 = get_J_0(orbital_operator_2e_prefactored_ss, eps, phi);
    const double J_1 = get_J_1(orbital_operator_2e_prefactored_ss, eps, phi);
    const double K = get_K(orbital_operator_2e_pure);
    const double K_0 = get_K_0(orbital_operator_2e_pure, eps, phi);
    const double K_1 = get_K_1(orbital_operator_2e_pure, eps, phi);
    const double L = get_L(orbital_operator_1e_pure);
    const double L_1 = get_L_1(orbital_operator_1e_pure, eps, phi);
    // Retutn:
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
