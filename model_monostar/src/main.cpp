#include <model_monostar/hardcoded_example.hpp>
#include <model_monostar/lin_alg.hpp>
#include <model_monostar/monostar_basis.hpp>
#include <model_monostar/monostar_hamiltonian_k0.hpp>
#include <model_monostar/monostar_kstate.hpp>
#include <model_monostar/monostar_site_state.hpp>

#include <extensions/adaptors.hpp>

#include <armadillo>

#include <iostream>
#include <iterator>

#include <cassert>  // TODO:  move to hamiltonian sourcefile
#include <complex>

using namespace std::complex_literals;

// #######################################################################
// ## DynamicMonostarHamiltonian                                        ##
// #######################################################################

namespace model_monostar {

class DynamicMonostarHamiltonian {
   public:
    DynamicMonostarHamiltonian(const size_t n_sites);
    void push_back_conjugated_states_to_basis(
        const DynamicMonostarUniqueKstate& generator,
        DynamicMonostarUniqueKstateBasis& basis) const;
    // ---- for: k_n =0
    void fill_k0_hamiltonian_matrix_coll(
        const DynamicMonostarUniqueKstateBasis& basis,
        size_t n_col,
        arma::mat& k0_hamiltonian_matrix) const;
    void fill_k0_hamiltonian_matrix(
        const DynamicMonostarUniqueKstateBasis& basis,
        arma::mat& k0_hamiltonian_matrix) const;
    arma::mat make_k0_hamiltonian_matrix(
        const DynamicMonostarUniqueKstateBasis& basis) const;
    // ---- for arbitrary k_n
    void fill_kn_hamiltonian_matrix_coll(
        const DynamicMonostarUniqueKstateBasis& basis,
        size_t n_col,
        arma::sp_cx_mat& kn_hamiltonian_matrix,
        const unsigned k_n) const;
    void fill_kn_hamiltonian_matrix(
        const DynamicMonostarUniqueKstateBasis& basis,
        arma::sp_cx_mat& kn_hamiltonian_matrix,
        const unsigned k_n) const;
    arma::sp_cx_mat make_kn_hamiltonian_matrix(
        const DynamicMonostarUniqueKstateBasis& basis,
        const unsigned k_n) const;

   private:
    const size_t _n_sites;
};

inline DynamicMonostarHamiltonian::DynamicMonostarHamiltonian(const size_t n_sites)
    : _n_sites(n_sites) {}

inline void DynamicMonostarHamiltonian::push_back_conjugated_states_to_basis(
    const DynamicMonostarUniqueKstate& generator,
    DynamicMonostarUniqueKstateBasis& basis) const {
    assert(generator.n_sites() == _n_sites);
    const auto generator_range = generator.to_range();
    for (size_t i = 0, j = 1; i < _n_sites; i++, j = (i + 1) % _n_sites) {
        if (*std::next(std::begin(generator_range), i) == gs &&
            *std::next(std::begin(generator_range), j) == gs) {
            const auto conjugated_range = generator_range | extension::boost::adaptors::refined(i, es) | extension::boost::adaptors::refined(j, es);
            const auto conjugated_kstate_ptr = std::make_shared<DynamicMonostarUniqueKstate>(conjugated_range, kstate::ctr_from_range);
            basis.add_element(conjugated_kstate_ptr);
        }
    }
    for (size_t i = 0, j = 1; i < _n_sites; i++, j = (i + 1) % _n_sites) {
        if (*std::next(std::begin(generator_range), i) == es &&
            *std::next(std::begin(generator_range), j) == es) {
            const auto conjugated_range = generator_range | extension::boost::adaptors::refined(i, gs) | extension::boost::adaptors::refined(j, gs);
            const auto conjugated_kstate_ptr = std::make_shared<DynamicMonostarUniqueKstate>(conjugated_range, kstate::ctr_from_range);
            basis.add_element(conjugated_kstate_ptr);
        }
    }
}

/// ------------------------------ for k_n = 0:

inline void DynamicMonostarHamiltonian::fill_k0_hamiltonian_matrix_coll(
    const DynamicMonostarUniqueKstateBasis& basis,
    const size_t ket_kstate_idx,
    arma::mat& k0_hamiltonian_matrix) const {
    assert(k0_hamiltonian_matrix.n_cols == basis.size());
    assert(k0_hamiltonian_matrix.n_rows == basis.size());
    const auto ket_kstate_ptr = basis.vec_index()[ket_kstate_idx];
    assert(ket_kstate_ptr);
    const auto& ket_kstate = (*ket_kstate_ptr).to_range();
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        if (*std::next(std::begin(ket_kstate), n_delta) == gs &&
            *std::next(std::begin(ket_kstate), n_delta_p1) == gs) {
            const auto bra_kstate = ket_kstate | extension::boost::adaptors::refined(n_delta, es) | extension::boost::adaptors::refined(n_delta_p1, es);
            if (const auto& bra_kstate_optional_idx = basis.find_element_and_get_its_ra_index(kstate::make_unique_shift(bra_kstate))) {
                const auto bra_kstate_idx = *bra_kstate_optional_idx;
                double pre_norm = _n_sites * (_n_sites / basis.vec_index()[bra_kstate_idx]->n_least_replication_shift()) * basis.vec_index()[bra_kstate_idx]->norm_factor() * basis.vec_index()[ket_kstate_idx]->norm_factor();
                k0_hamiltonian_matrix(bra_kstate_idx, ket_kstate_idx) = k0_hamiltonian_matrix(bra_kstate_idx, ket_kstate_idx) + pre_norm * 0.5;
            }
        }
    }
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        if (*std::next(std::begin(ket_kstate), n_delta) == es &&
            *std::next(std::begin(ket_kstate), n_delta_p1) == es) {
            const auto bra_kstate = ket_kstate | extension::boost::adaptors::refined(n_delta, gs) | extension::boost::adaptors::refined(n_delta_p1, gs);
            if (const auto& bra_kstate_optional_idx = basis.find_element_and_get_its_ra_index(kstate::make_unique_shift(bra_kstate))) {
                const auto bra_kstate_idx = *bra_kstate_optional_idx;
                double pre_norm = _n_sites * (_n_sites / basis.vec_index()[bra_kstate_idx]->n_least_replication_shift()) * basis.vec_index()[bra_kstate_idx]->norm_factor() * basis.vec_index()[ket_kstate_idx]->norm_factor();
                k0_hamiltonian_matrix(bra_kstate_idx, ket_kstate_idx) = k0_hamiltonian_matrix(bra_kstate_idx, ket_kstate_idx) + pre_norm * 0.5;
            }
        }
    }
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        const bool is_the_same = *std::next(std::begin(ket_kstate), n_delta) == *std::next(std::begin(ket_kstate), n_delta_p1);
        double pre_norm = _n_sites * (_n_sites / basis.vec_index()[ket_kstate_idx]->n_least_replication_shift()) * basis.vec_index()[ket_kstate_idx]->norm_factor() * basis.vec_index()[ket_kstate_idx]->norm_factor();
        k0_hamiltonian_matrix(ket_kstate_idx, ket_kstate_idx) += pre_norm * (is_the_same ? -1.0 / 4.0 : +1.0 / 4.0);
    }
}

inline void DynamicMonostarHamiltonian::fill_k0_hamiltonian_matrix(
    const DynamicMonostarUniqueKstateBasis& basis,
    arma::mat& k0_hamiltonian_matrix) const {
    assert(k0_hamiltonian_matrix.n_cols == basis.size());
    assert(k0_hamiltonian_matrix.n_rows == basis.size());
    for (arma::uword ket_kstate_idx = 0; ket_kstate_idx < basis.size(); ket_kstate_idx++) {
        fill_k0_hamiltonian_matrix_coll(basis, ket_kstate_idx, k0_hamiltonian_matrix);
    }
}

inline arma::mat
DynamicMonostarHamiltonian::make_k0_hamiltonian_matrix(
    const DynamicMonostarUniqueKstateBasis& basis) const {
    arma::mat k0_hamiltonian_matrix(basis.size(), basis.size(), arma::fill::zeros);
    fill_k0_hamiltonian_matrix(basis, k0_hamiltonian_matrix);
    return k0_hamiltonian_matrix;
}

/// ------------------------------ arbitrary k_n:

inline void DynamicMonostarHamiltonian::fill_kn_hamiltonian_matrix_coll(
    const DynamicMonostarUniqueKstateBasis& basis,
    const size_t ket_kstate_idx,
    arma::sp_cx_mat& kn_hamiltonian_matrix,
    const unsigned k_n) const {
    assert(kn_hamiltonian_matrix.n_cols == basis.size());
    assert(kn_hamiltonian_matrix.n_rows == basis.size());
    const auto ket_kstate_ptr = basis.vec_index()[ket_kstate_idx];
    assert(ket_kstate_ptr);
    const auto& ket_kstate = (*ket_kstate_ptr).to_range();
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        if (*std::next(std::begin(ket_kstate), n_delta) == gs &&
            *std::next(std::begin(ket_kstate), n_delta_p1) == gs) {
            const auto bra_kstate = ket_kstate | extension::boost::adaptors::refined(n_delta, es) | extension::boost::adaptors::refined(n_delta_p1, es);
            if (const auto& bra_kstate_optional_idx = basis.find_element_and_get_its_ra_index(kstate::make_unique_shift(bra_kstate))) {
                const auto bra_kstate_idx = *bra_kstate_optional_idx;
                double pre_norm = _n_sites * basis.vec_index()[bra_kstate_idx]->norm_factor() * basis.vec_index()[ket_kstate_idx]->norm_factor();
                std::complex<double> sum_phase_factors = 0;
                const size_t bra_n_unique_shift = kstate::n_unique_shift(bra_kstate);
                const size_t bra_n_least_replication_shift = basis.vec_index()[bra_kstate_idx]->n_least_replication_shift();
                const size_t bra_n_replicas = _n_sites / bra_n_least_replication_shift;
                for (unsigned n_contribution = 0; n_contribution < bra_n_replicas; n_contribution++) {
                    // sign in front of bra_n_unique_shift in the below formula has to be check!
                    const double exponent = 2 * arma::datum::pi * k_n / _n_sites * (-(int)bra_n_unique_shift + (int)n_contribution * (int)bra_n_least_replication_shift);
                    sum_phase_factors += std::exp(1.0i * exponent);
                }
                kn_hamiltonian_matrix(bra_kstate_idx, ket_kstate_idx) += pre_norm * sum_phase_factors * 0.5;
            }
        }
    }
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        if (*std::next(std::begin(ket_kstate), n_delta) == es &&
            *std::next(std::begin(ket_kstate), n_delta_p1) == es) {
            const auto bra_kstate = ket_kstate | extension::boost::adaptors::refined(n_delta, gs) | extension::boost::adaptors::refined(n_delta_p1, gs);
            if (const auto& bra_kstate_optional_idx = basis.find_element_and_get_its_ra_index(kstate::make_unique_shift(bra_kstate))) {
                const auto bra_kstate_idx = *bra_kstate_optional_idx;
                double pre_norm = _n_sites * basis.vec_index()[bra_kstate_idx]->norm_factor() * basis.vec_index()[ket_kstate_idx]->norm_factor();
                std::complex<double> sum_phase_factors = 0;
                const size_t bra_n_unique_shift = kstate::n_unique_shift(bra_kstate);
                const size_t bra_n_least_replication_shift = basis.vec_index()[bra_kstate_idx]->n_least_replication_shift();
                const size_t bra_n_replicas = _n_sites / bra_n_least_replication_shift;
                for (unsigned n_contribution = 0; n_contribution < bra_n_replicas; n_contribution++) {
                    const double exponent = 2 * arma::datum::pi * k_n / _n_sites * (-(int)bra_n_unique_shift + (int)n_contribution * (int)bra_n_least_replication_shift);
                    sum_phase_factors += std::exp(1.0i * exponent);
                    std::cout << "IDXs " << bra_kstate_idx << " x " << ket_kstate_idx << std::endl;
                    std::cout << bra_kstate << " x " << ket_kstate << std::endl;
                    std::cout << "n_contribution " << n_contribution << ", "
                              << "bra_n_unique_shift " << bra_n_unique_shift << ", "
                              << "exponent/pi " << 2.0 * k_n / _n_sites * (-(int)bra_n_unique_shift + (int)n_contribution * (int)bra_n_least_replication_shift) << ", "
                              << std::endl;
                }
                kn_hamiltonian_matrix(bra_kstate_idx, ket_kstate_idx) += pre_norm * sum_phase_factors * 0.5;
            }
        }
    }
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        const bool is_the_same = *std::next(std::begin(ket_kstate), n_delta) == *std::next(std::begin(ket_kstate), n_delta_p1);
        double pre_norm = _n_sites * (_n_sites / basis.vec_index()[ket_kstate_idx]->n_least_replication_shift()) * basis.vec_index()[ket_kstate_idx]->norm_factor() * basis.vec_index()[ket_kstate_idx]->norm_factor();
        kn_hamiltonian_matrix(ket_kstate_idx, ket_kstate_idx) += pre_norm * (is_the_same ? -1.0 / 4.0 : +1.0 / 4.0);
    }
}

inline void DynamicMonostarHamiltonian::fill_kn_hamiltonian_matrix(
    const DynamicMonostarUniqueKstateBasis& basis,
    arma::sp_cx_mat& kn_hamiltonian_matrix,
    const unsigned k_n) const {
    assert(kn_hamiltonian_matrix.n_cols == basis.size());
    assert(kn_hamiltonian_matrix.n_rows == basis.size());
    for (arma::uword ket_kstate_idx = 0; ket_kstate_idx < basis.size(); ket_kstate_idx++) {
        fill_kn_hamiltonian_matrix_coll(basis, ket_kstate_idx, kn_hamiltonian_matrix, k_n);
    }
}

inline arma::sp_cx_mat
DynamicMonostarHamiltonian::make_kn_hamiltonian_matrix(
    const DynamicMonostarUniqueKstateBasis& basis,
    const unsigned k_n) const {
    arma::sp_cx_mat kn_hamiltonian_matrix(basis.size(), basis.size());
    fill_kn_hamiltonian_matrix(basis, kn_hamiltonian_matrix, k_n);
    return kn_hamiltonian_matrix;
}

}  // namespace model_monostar

// #######################################################################
// ## main...                                                           ##
// #######################################################################

void populate_pt_basis(
    const model_monostar::DynamicMonostarHamiltonian& hamiltonian,
    const unsigned max_pt_order,
    model_monostar::DynamicMonostarUniqueKstateBasis& basis) {
    // assert();
    unsigned last_chunk_size = basis.size();
    for (unsigned pt_order = 0; pt_order < max_pt_order; pt_order++) {
        const unsigned old_size = basis.size();
        assert(last_chunk_size <= old_size);
        for (unsigned idx = old_size - last_chunk_size; idx < old_size; idx++) {
            const auto el = basis.vec_index()[idx];
            assert(el);
            hamiltonian.push_back_conjugated_states_to_basis(*el, basis);
        }
        const unsigned new_size = basis.size();
        last_chunk_size = new_size - old_size;
    }
}

double bfpt_gs(const size_t n_sites, const unsigned max_pt_order) {
    // Define hamiltonian and basis:
    model_monostar::DynamicMonostarUniqueKstateBasis basis{n_sites};
    model_monostar::DynamicMonostarHamiltonian hamiltonian{n_sites};
    // Define pt_orders subspace basis:
    basis.add_element(std::make_shared<model_monostar::DynamicMonostarUniqueKstate>(model_monostar::classical_gs_kstate(n_sites)));
    // ----
    std::cout << basis;
    // ----
    // Generate higher pt_orders subspace basis:
    populate_pt_basis(hamiltonian, max_pt_order, basis);
    // ----
    std::cout << basis;
    // ----
    const auto k0_hamiltonian_matrix = hamiltonian.make_k0_hamiltonian_matrix(basis);
    // ----
    // std::cout << k0_hamiltonian_matrix;
    // ----
    arma::vec eigen_values;
    arma::mat eigen_vectors;
    arma::eig_sym(eigen_values, eigen_vectors, k0_hamiltonian_matrix);
    // ----
    std::cout << eigen_values;
    std::cout << "min eigen_value: " << eigen_values(0) << std::endl;
    // ----
    return eigen_values(0);
}

double bfpt_goldston(const size_t n_sites, const unsigned max_pt_order) {
    // Define hamiltonian and basis:
    model_monostar::DynamicMonostarUniqueKstateBasis basis{n_sites};
    model_monostar::DynamicMonostarHamiltonian hamiltonian{n_sites};
    // Define pt_orders subspace basis:
    basis.add_element(std::make_shared<model_monostar::DynamicMonostarUniqueKstate>(model_monostar::classical_es_kstate(n_sites)));
    // ----
    std::cout << basis;
    // ----
    // Generate higher pt_orders subspace basis:
    populate_pt_basis(hamiltonian, max_pt_order, basis);
    // ----
    std::cout << basis;
    // ----
    const auto k0_hamiltonian_matrix = hamiltonian.make_k0_hamiltonian_matrix(basis);
    // ----
    // std::cout << k0_hamiltonian_matrix;
    // ----
    arma::vec eigen_values;
    arma::mat eigen_vectors;
    arma::eig_sym(eigen_values, eigen_vectors, k0_hamiltonian_matrix);
    // ----
    std::cout << eigen_values;
    std::cout << "min eigen_value: " << eigen_values(0) << std::endl;
    // ----
    return eigen_values(0);
}

double bfpt_kn_es(const size_t n_sites, const unsigned max_pt_order, const unsigned k_n) {
    // Define hamiltonian and basis:
    model_monostar::DynamicMonostarUniqueKstateBasis basis{n_sites};
    model_monostar::DynamicMonostarHamiltonian hamiltonian{n_sites};
    // Define pt_orders subspace basis:
    basis.add_element(std::make_shared<model_monostar::DynamicMonostarUniqueKstate>(model_monostar::classical_es_kstate(n_sites)));
    // ----
    std::cout << basis;
    // ----
    // Generate higher pt_orders subspace basis:
    populate_pt_basis(hamiltonian, max_pt_order, basis);
    // ----
    std::cout << basis;
    // ----
    const auto kn_hamiltonian_matrix = hamiltonian.make_kn_hamiltonian_matrix(basis, k_n);
    // ----
    // std::cout << kn_hamiltonian_matrix;
    // ----
    arma::vec eigen_values;
    arma::cx_mat eigen_vectors;
    arma::eig_sym(eigen_values, eigen_vectors, arma::cx_mat(kn_hamiltonian_matrix));
    // ----
    std::cout << "**** arma lin_alg:" << std::endl;
    std::cout << eigen_values;
    std::cout << "min eigen_value: " << eigen_values(0) << std::endl;
    std::cout << "eigen_vectors.col(0): " << std::endl
              << eigen_vectors.col(0) / eigen_vectors(0, 0) << std::endl;
    // std::cout << "eigen_vectors.col(1): " << std::endl
    //           << eigen_vectors.col(1) / eigen_vectors(0, 1) << std::endl;
    // std::cout << "eigen_vectors.col(5): " << std::endl
    //           << eigen_vectors.col(5) / eigen_vectors(0, 5) << std::endl;
    // ----
    std::cout << "**** my lin_alg spike:" << std::endl;
    arma::vec eigen_values_2;
    arma::cx_mat eigen_vectors_2;
    lin_alg::eigs_sym(eigen_values_2, eigen_vectors_2, kn_hamiltonian_matrix, 1, "sa", 1e-6);
    std::cout << eigen_values_2;
    std::cout << "min eigen_value: " << eigen_values_2(0) << std::endl;
    std::cout << "eigen_vectors.col(0): " << std::endl
              << eigen_vectors_2.col(0) / eigen_vectors_2(0, 0) << std::endl;
    // std::cout << "eigen_vectors.col(1): " << std::endl
    //           << eigen_vectors_2.col(1) / eigen_vectors_2(0, 1) << std::endl;
    // std::cout << "eigen_vectors.col(5): " << std::endl
    //           << eigen_vectors_2.col(5) / eigen_vectors_2(0, 5) << std::endl;
    // ----
    return eigen_values_2(0);
}

int main() {
    const unsigned max_pt_order = 1;
    const size_t n_sites = 8;
    const unsigned k_n = 1;
    std::cout << "------------------------------------------" << std::endl;
    double gs_energy = bfpt_gs(n_sites, max_pt_order);
    std::cout << "------------------------------------------" << std::endl;
    double es_gold_energy = bfpt_goldston(n_sites, max_pt_order);
    std::cout << "------------------------------------------" << std::endl;
    double es_energy = bfpt_kn_es(n_sites, max_pt_order, k_n);
    std::cout << "------------------------------------------" << std::endl;
    std::cout << " gs enery = " << gs_energy << std::endl;
    std::cout << " es goldston = " << es_gold_energy << std::endl;
    std::cout << "    excitation enery goldston = " << es_gold_energy - gs_energy << std::endl;
    std::cout << " k_n = " << k_n << std::endl;
    std::cout << " es enery (k_n) = " << es_energy << std::endl;
    std::cout << "    excitation enery (k_n) = " << es_energy - gs_energy << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    //model_monostar::do_hardcoded_example_analyse();
    return 0;
}
