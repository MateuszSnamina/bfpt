#include <model_monostar/hardcoded_example.hpp>
#include <model_monostar/monostar_basis.hpp>
#include <model_monostar/monostar_hamiltonian_k0.hpp>
#include <model_monostar/monostar_kstate.hpp>
#include <model_monostar/monostar_site_state.hpp>

#include <linear_algebra/linear_algebra.hpp>

#include <kstate/basis_populate.hpp>

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

class DynamicMonostarHamiltonian : public kstate::IDynamicUniqueKstatePopulator<MonostarSiteState> {
   public:
    DynamicMonostarHamiltonian(const size_t n_sites);
    // Generates all conjugated states:
    void push_back_coupled_states_to_basis(
        const DynamicMonostarUniqueKstate& generator,
        DynamicMonostarUniqueKstateBasis& basis) const override;
    // Generates hamiltonian matrix:
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

inline void DynamicMonostarHamiltonian::push_back_coupled_states_to_basis(
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

inline void DynamicMonostarHamiltonian::fill_kn_hamiltonian_matrix_coll(
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
        if (*std::next(std::begin(ket_kstate), n_delta) == gs &&
            *std::next(std::begin(ket_kstate), n_delta_p1) == gs) {
            const auto bra_kstate =
                ket_kstate | extension::boost::adaptors::refined(n_delta, es) | extension::boost::adaptors::refined(n_delta_p1, es);
            if (const auto& bra_kstate_optional_idx = basis.find_element_and_get_its_ra_index(kstate::make_unique_shift(bra_kstate))) {
                const auto bra_kstate_idx = *bra_kstate_optional_idx;
                double pre_norm = _n_sites * basis.vec_index()[bra_kstate_idx]->norm_factor() * basis.vec_index()[ket_kstate_idx]->norm_factor();
                std::complex<double> sum_phase_factors = 0;
                const size_t bra_n_unique_shift = kstate::n_unique_shift(bra_kstate);
                const size_t bra_n_least_replication_shift = basis.vec_index()[bra_kstate_idx]->n_least_replication_shift();
                const size_t bra_n_replicas = _n_sites / bra_n_least_replication_shift;
                for (unsigned n_contribution = 0; n_contribution < bra_n_replicas; n_contribution++) {
                    // sign in front of bra_n_unique_shift in the below formula has to be check!
                    const int exponent_n = -(int)bra_n_unique_shift + (int)n_contribution * (int)bra_n_least_replication_shift;
                    const double exponent_r = 2 * arma::datum::pi * k_n / _n_sites * exponent_n;
                    sum_phase_factors += std::exp(1.0i * exponent_r);
                }
                kn_hamiltonian_matrix(bra_kstate_idx, ket_kstate_idx) += pre_norm * sum_phase_factors * 0.5;
                kn_hamiltonian_matrix(ket_kstate_idx, bra_kstate_idx) += std::conj(pre_norm * sum_phase_factors * 0.5);
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

struct CommonRecipePrintFlags {
    bool print_unpopulated_basis_flag = false;
    bool print_unpopulated_basis_size_flag = true;
    bool print_populated_basis_flag = false;
    bool print_populated_basis_size_flag = true;
    bool print_sp_hamiltonian_flag = false;
    bool print_hamiltonian_flag = false;
    bool print_eigen_values_flag = true;
    bool print_eigen_vectors_flag = false;
};

double
do_common_recipe(model_monostar::DynamicMonostarUniqueKstateBasis& basis, const unsigned max_pt_order, const unsigned k_n,
                 CommonRecipePrintFlags print_flags) {
    arma::wall_clock timer;
    const size_t n_sites = basis.n_sites();
    const std::string message_prefix = "[common-recipe] ";
    const std::string progress_tag = "[progress] ";
    const std::string data_tag = "[data    ] ";
    const std::string time_tag = "[time    ] ";
    assert(k_n < n_sites);
    // --------------------------------------------------
    if (print_flags.print_unpopulated_basis_flag) {
        std::cout << message_prefix << data_tag << "Unpopulated basis (0'th pt-order basis):";
        std::cout << basis;
    }
    if (print_flags.print_unpopulated_basis_size_flag) {
        std::cout << message_prefix << data_tag << "Unpopulated basis (0'th pt-order basis) size: "
                  << basis.size() << "." << std::endl;
    }
    // --------------------------------------------------
    // Define hamiltonian and basis:
    std::cout << message_prefix << progress_tag << "About to populate pt-basis." << std::endl;
    model_monostar::DynamicMonostarHamiltonian hamiltonian{n_sites};
    // Generate higher pt-orders subspace basis:
    timer.tic();
    kstate::populate_pt_basis(hamiltonian, max_pt_order, basis);
    const double time_populating_pt_basis = timer.toc();
    std::cout << message_prefix << time_tag << "Populating pt-basis took " << time_populating_pt_basis << "s." << std::endl;
    std::cout << message_prefix << progress_tag << "Has populated pt-basis." << std::endl;
    // --------------------------------------------------
    if (print_flags.print_populated_basis_flag) {
        std::cout << message_prefix << data_tag << "Populated basis (" << max_pt_order << "'th pt-order basis):";
        std::cout << basis;
    }
    if (print_flags.print_populated_basis_size_flag) {
        std::cout << message_prefix << data_tag << "Populated basis (" << max_pt_order << "'th pt-order basis) size: "
                  << basis.size() << "." << std::endl;
    }
    // --------------------------------------------------
    // Generate hamiltoniam matrix:
    std::cout << message_prefix << progress_tag << "About to generate hamiltoniam." << std::endl;
    timer.tic();
    const auto kn_hamiltonian_matrix = hamiltonian.make_kn_hamiltonian_matrix(basis, k_n);
    const double time_generating_kn_hamiltonian_matrix = timer.toc();
    std::cout << message_prefix << time_tag << "Generating kn-hamiltoniam matrix took " << time_generating_kn_hamiltonian_matrix << "s." << std::endl;
    std::cout << message_prefix << progress_tag << "Has generated kn-hamiltoniam." << std::endl;
    // --------------------------------------------------
    if (print_flags.print_sp_hamiltonian_flag) {
        std::cout << message_prefix << data_tag << "kn_hamiltonian_matrix:";
        std::cout << kn_hamiltonian_matrix;
    }
    if (print_flags.print_hamiltonian_flag) {
        std::cout << message_prefix << data_tag << "kn_hamiltonian_matrix:";
        std::cout << arma::cx_mat(kn_hamiltonian_matrix);
    }
    // --------------------------------------------------
    // arma::vec eigen_values_debug;
    // arma::cx_mat eigen_vectors_debug;
    // arma::eig_sym(eigen_values_debug, eigen_vectors_debug, arma::cx_mat(kn_hamiltonian_matrix));
    // std::cout << eigen_values_debug;
    // std::cout << "min eigen_value: " << eigen_values_debug(0) << std::endl;
    // std::cout << "eigen_vectors_debug.col(0): " << std::endl
    //           << eigen_vectors_debug.col(0) / eigen_vectors_debug(0, 0) << std::endl;
    // --------------------------------------------------
    std::cout << message_prefix << progress_tag << "About to solve eigen problem." << std::endl;
    timer.tic();
    const auto& eigs_sym_result = lin_alg::eigs_sym(kn_hamiltonian_matrix, 1, 3, "sa", 1e-6);
    const double time_solving_eigen_problem = timer.toc();
    if (eigs_sym_result.is_err()) {
        std::cout << message_prefix << time_tag << "Solving eigen problem took: " << time_solving_eigen_problem << "s." << std::endl;
        std::cout << message_prefix << progress_tag << "Failed to solve eigen problem." << std::endl;
        return arma::datum::nan;
    }
    const auto eigen_info = eigs_sym_result.unwrap();
    // const auto eigen_info = lin_alg::eigs_sym(kn_hamiltonian_matrix, 1, 3, "sa", 1e-6).unwrap();
    const arma::vec& eigen_values = eigen_info.eigen_values;
    const arma::cx_mat& eigen_vectors = eigen_info.eigen_vectors;
    std::cout << message_prefix << time_tag << "Solving eigen problem took: " << time_solving_eigen_problem << "s." << std::endl;
    std::cout << message_prefix << progress_tag << "Has solved eigen problem." << std::endl;
    // --------------------------------------------------
    if (print_flags.print_eigen_values_flag) {
        std::cout << message_prefix << data_tag << "eigen_values:" << std::endl;
        std::cout << message_prefix << eigen_values;
        // std::cout << "min eigen_value: " << eigen_values(0) << std::endl;
    }
    if (print_flags.print_eigen_vectors_flag) {
        std::cout << message_prefix << data_tag << "eigen_vectors:" << std::endl;
        std::cout << message_prefix << eigen_vectors;
        // std::cout << "eigen_vectors.col(0): " << std::endl
        //           << eigen_vectors.col(0) / eigen_vectors(0, 0) << std::endl;
    }
    // --------------------------------------------------
    return eigen_values(0);
}

double bfpt_gs(const size_t n_sites, const unsigned max_pt_order) {
    CommonRecipePrintFlags print_flags;
    model_monostar::DynamicMonostarUniqueKstateBasis basis{n_sites};
    basis.add_element(std::make_shared<model_monostar::DynamicMonostarUniqueKstate>(model_monostar::classical_gs_kstate(n_sites)));
    return do_common_recipe(basis, max_pt_order, 0, print_flags);
}

double bfpt_kn_es(const size_t n_sites, const unsigned max_pt_order, const unsigned k_n) {
    CommonRecipePrintFlags print_flags;
    model_monostar::DynamicMonostarUniqueKstateBasis basis{n_sites};
    basis.add_element(std::make_shared<model_monostar::DynamicMonostarUniqueKstate>(model_monostar::classical_es_kstate(n_sites)));
    return do_common_recipe(basis, max_pt_order, k_n, print_flags);
}

double bfpt_goldston(const size_t n_sites, const unsigned max_pt_order) {
    return bfpt_kn_es(n_sites, max_pt_order, 0);
}

int main() {
    // const unsigned max_pt_order = 1;
    // const size_t n_sites = 8;
    // const unsigned k_n = 1;
    const unsigned max_pt_order = 6;
    const size_t n_sites = 20;
    const unsigned k_n = 1;
    std::cout << "------------------------------------------" << std::endl;
    double gs_energy = bfpt_gs(n_sites, max_pt_order);
    std::cout << "------------------------------------------" << std::endl;
    double es_gold_energy = bfpt_goldston(n_sites, max_pt_order);
    std::cout << "------------------------------------------" << std::endl;
    double es_energy = bfpt_kn_es(n_sites, max_pt_order, k_n);
    std::cout << "------------------------------------------" << std::endl;
    std::cout << " max_pt_order = " << max_pt_order << std::endl;
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
