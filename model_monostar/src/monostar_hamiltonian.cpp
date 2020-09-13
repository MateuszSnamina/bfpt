#include <bfpt_common/do_common_recipie.hpp>  // TODO in future move down.
#include <bfpt_common/populate_pt_basis.hpp>  // TODO in future move down.

#include <model_monostar/monostar_hamiltonian.hpp>

#include <extensions/adaptors.hpp>

#include <cassert>
#include <complex>

using namespace std::complex_literals;

// #######################################################################
// ## DynamicMonostarHamiltonian                                        ##
// #######################################################################

namespace model_monostar {

DynamicMonostarHamiltonian::DynamicMonostarHamiltonian(const size_t n_sites)
    : _n_sites(n_sites) {}

void DynamicMonostarHamiltonian::push_back_coupled_states_to_basis(
        const DynamicMonostarUniqueKstate& generator,
        DynamicMonostarUniqueKstateBasis& basis) const {
    assert(generator.n_sites() == _n_sites);
    const auto generator_range = generator.to_range();
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        if (*std::next(std::begin(generator_range), n_delta) == gs &&
                *std::next(std::begin(generator_range), n_delta_p1) == gs) {
            const auto conjugated_range =
                    generator_range | extension::boost::adaptors::refined(n_delta, es) | extension::boost::adaptors::refined(n_delta_p1, es);
            const auto conjugated_kstate_ptr = std::make_shared<DynamicMonostarUniqueKstate>(conjugated_range, kstate::ctr_from_range);
            basis.add_element(conjugated_kstate_ptr);
        }
    }
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        if (*std::next(std::begin(generator_range), n_delta) == es &&
                *std::next(std::begin(generator_range), n_delta_p1) == es) {
            const auto conjugated_range =
                    generator_range | extension::boost::adaptors::refined(n_delta, gs) | extension::boost::adaptors::refined(n_delta_p1, gs);
            const auto conjugated_kstate_ptr = std::make_shared<DynamicMonostarUniqueKstate>(conjugated_range, kstate::ctr_from_range);
            basis.add_element(conjugated_kstate_ptr);
        }
    }
}

/// ------------------------------ arbitrary k_n:

void DynamicMonostarHamiltonian::fill_kn_hamiltonian_matrix_coll(
        const DynamicMonostarUniqueKstateBasis& basis,
        const size_t ket_kstate_idx,
        arma::sp_cx_mat& kn_hamiltonian_matrix,
        const unsigned k_n) const {
    assert(kn_hamiltonian_matrix.n_cols == kn_hamiltonian_matrix.n_rows);
    assert(kn_hamiltonian_matrix.n_rows == basis.size());
    const auto ket_kstate_ptr = basis.vec_index()[ket_kstate_idx];
    assert(ket_kstate_ptr);
    const auto& ket_kstate = (*ket_kstate_ptr).to_range();
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        const auto ket_site_1 = *std::next(std::begin(ket_kstate), n_delta);
        const auto ket_site_2 = *std::next(std::begin(ket_kstate), n_delta_p1);
        if (ket_site_1 == gs && ket_site_2 == gs) {
            const auto bra_site_1 = es;
            const auto bra_site_2 = es;
            const auto bra_kstate = ket_kstate
                    | extension::boost::adaptors::refined(n_delta, bra_site_1)
                    | extension::boost::adaptors::refined(n_delta_p1, bra_site_2);
            if (const auto& bra_kstate_optional_idx = basis.find_element_and_get_its_ra_index(kstate::make_unique_shift(bra_kstate))) {
                const auto bra_kstate_idx = *bra_kstate_optional_idx;
                double pre_norm_1 = _n_sites * basis.vec_index()[bra_kstate_idx]->norm_factor() * basis.vec_index()[ket_kstate_idx]->norm_factor();
                const size_t bra_n_unique_shift = kstate::n_unique_shift(bra_kstate);
                const size_t bra_n_least_replication_shift = basis.vec_index()[bra_kstate_idx]->n_least_replication_shift();
                const size_t bra_n_replicas = _n_sites / bra_n_least_replication_shift;
                const int exponent_n = (int)bra_n_unique_shift;
                const double exponent_r = 2 * arma::datum::pi * k_n / _n_sites * exponent_n;
                const std::complex<double> neo_sum_phase_factors =
                        std::exp(1.0i * exponent_r) * (bra_n_replicas == 1 ? 1.0 : (k_n == 0 ? bra_n_replicas : 0.0));
                const std::complex<double> pre_norm_2 = pre_norm_1 * neo_sum_phase_factors;
                kn_hamiltonian_matrix(bra_kstate_idx, ket_kstate_idx) += pre_norm_2 * 0.5;
                kn_hamiltonian_matrix(ket_kstate_idx, bra_kstate_idx) += std::conj(pre_norm_2 * 0.5);
            }
        }
    }
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        const auto ket_site_1 = *std::next(std::begin(ket_kstate), n_delta);
        const auto ket_site_2 = *std::next(std::begin(ket_kstate), n_delta_p1);
        const bool is_the_same = ket_site_2 == ket_site_1;
        const double pre_norm_1 = _n_sites * ket_kstate_ptr->norm_factor() * ket_kstate_ptr->norm_factor();
        const double pre_norm_2 = (_n_sites / ket_kstate_ptr->n_least_replication_shift()) * pre_norm_1;
        kn_hamiltonian_matrix(ket_kstate_idx, ket_kstate_idx) += pre_norm_2 * (is_the_same ? -1.0 / 4.0 : +1.0 / 4.0);
    }
}

void DynamicMonostarHamiltonian::fill_kn_hamiltonian_matrix(
        const DynamicMonostarUniqueKstateBasis& basis,
        arma::sp_cx_mat& kn_hamiltonian_matrix,
        const unsigned k_n) const {
    assert(kn_hamiltonian_matrix.n_cols == basis.size());
    assert(kn_hamiltonian_matrix.n_rows == basis.size());
    for (arma::uword ket_kstate_idx = 0; ket_kstate_idx < basis.size(); ket_kstate_idx++) {
        fill_kn_hamiltonian_matrix_coll(basis, ket_kstate_idx, kn_hamiltonian_matrix, k_n);
    }
}

arma::sp_cx_mat
DynamicMonostarHamiltonian::make_kn_hamiltonian_matrix(
        const DynamicMonostarUniqueKstateBasis& basis,
        const unsigned k_n) const {
    arma::sp_cx_mat kn_hamiltonian_matrix(basis.size(), basis.size());
    fill_kn_hamiltonian_matrix(basis, kn_hamiltonian_matrix, k_n);
    return kn_hamiltonian_matrix;
}

}  // namespace model_monostar
