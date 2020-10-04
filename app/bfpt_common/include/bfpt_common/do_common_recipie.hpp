#ifndef BFPT_COMMON_DO_COMMON_RECIPIE_HPP
#define BFPT_COMMON_DO_COMMON_RECIPIE_HPP

#include <bfpt_common/i_kstate_operator_matrix.hpp>
#include <bfpt_common/i_kstate_basis_populator.hpp>
#include <bfpt_common/generate_populated_basis.hpp> // TODO change to  <bfpt_common/generate_basis.hpp>
#include <bfpt_common/generate_operator_matrix.hpp>
#include <bfpt_common/common_recipe_print_flags.hpp>
#include <bfpt_common/calculate_reduced_density_operator.hpp>
#include <bfpt_common/density_operator.hpp>

#include <linear_algebra/linear_algebra.hpp>

#include <kstate/basis.hpp>
#include <kstate/kstate_concrete.hpp>

#include <utility/result.hpp>

#include <extensions/stream_fromat_stacker.hpp>

#include <armadillo>

#include <iostream>
#include <string>

#include <cassert>
#include <iomanip>
#include <algorithm>
#include <complex>
#include <cmath>

namespace bfpt_common {

const std::string message_prefix = "[common-recipe] ";
const std::string progress_tag = "[progress] ";
const std::string data_tag = "[data    ] ";
const std::string time_tag = "[time    ] ";

}

// #######################################################################
// ## print_density_operator                                            ##
// #######################################################################

namespace bfpt_common {

template<typename SiteStateT>
void pretty_density_operator_12(
        const DensityOperator12<SiteStateT>& density_operator,
        std::string print_outer_prefix = "") {
    // Stream RAII:
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    // Print:
    std::cout << print_outer_prefix << message_prefix << data_tag
              << "density operator 12:" << std::endl;
    for (const auto& _ : density_operator) {
        const std::pair<StateKernel12<SiteStateT>, StateKernel12<SiteStateT>> density_matrix_indices = _.first;
        const std::complex<double> value = _.second;
        const StateKernel12<SiteStateT> bra_kenrel = density_matrix_indices.first;
        const StateKernel12<SiteStateT> ket_kenrel = density_matrix_indices.second;
        const auto bra_site_0 = bra_kenrel.state_1;
        const auto bra_site_1 = bra_kenrel.state_2;
        const auto ket_site_0 = ket_kenrel.state_1;
        const auto ket_site_1 = ket_kenrel.state_2;
        std::cout << print_outer_prefix << message_prefix << data_tag
                  << "(" << bra_site_0 << "⊗" << bra_site_1 << ")" << " "
                  << "(" << ket_site_0 << "⊗" << ket_site_1 << ")" << " "
                  << std::showpos << value << std::endl;
    }
}

}

// #######################################################################
// ## pretty_print                                                      ##
// #######################################################################

namespace bfpt_common {

template<typename KstateT>
void pretty_print(
        kstate::Basis<KstateT>& basis,
        const arma::vec& eigen_values,
        const arma::cx_mat& eigen_vectors,
        std::pair<unsigned, unsigned> print_pretty_min_max_n_kstates,
        double print_pretty_probability_treshold,
        std::string print_outer_prefix = "") {
    // [[expect]]:
    assert(eigen_vectors.n_cols == eigen_values.n_rows);
    assert(basis.size() == eigen_vectors.n_rows);
    // helpers:
    const auto basis_size = basis.size();
    struct KstateIdxAndProbability {
        unsigned idx;
        double probability;
    };
    const auto order_kstateIdx_and_probability = []
            (const KstateIdxAndProbability& lhs, const KstateIdxAndProbability& rhs) -> bool {
        return lhs.probability > rhs.probability;
    };
    // Stream RAII:
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    // Print:
    for (unsigned n_eigien_vector = 0; n_eigien_vector < eigen_vectors.n_cols; n_eigien_vector++) {
        const auto eigen_value = eigen_values(n_eigien_vector);
        const auto eigien_vector = eigen_vectors.col(n_eigien_vector);
        std::cout << print_outer_prefix << message_prefix << data_tag
                  << "eigen vector no: " << std::setw(6) << n_eigien_vector << ", "
                  << "eigen energy: " << eigen_value
                  << "." << std::endl;
        const std::vector<KstateIdxAndProbability> kstateIdx_and_probability_vector =
                [&basis_size, &eigien_vector, &order_kstateIdx_and_probability]() {
            std::vector<KstateIdxAndProbability> kstateIdx_and_probability_vector_builder;
            kstateIdx_and_probability_vector_builder.resize(basis_size);
            for(unsigned idx = 0; idx < basis_size; idx++) {
                kstateIdx_and_probability_vector_builder[idx] = {idx, std::norm(eigien_vector(idx))};
            }
            std::stable_sort(std::begin(kstateIdx_and_probability_vector_builder), std::end(kstateIdx_and_probability_vector_builder), order_kstateIdx_and_probability);
            return kstateIdx_and_probability_vector_builder;
        }();
        const std::complex<double> presentation_layer_amplitude_factor = [&kstateIdx_and_probability_vector, &eigien_vector] () {
            const unsigned print_idx = 0;
            const unsigned idx = kstateIdx_and_probability_vector[print_idx].idx;
            const std::complex<double> amplitude = eigien_vector(idx);
            return 1.0 / (amplitude / std::abs(amplitude));
        }();
        const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
        std::cout << std::fixed << std::setprecision(9);
        double accumulated_probability = 0.0;
        for(unsigned print_idx = 0; print_idx < basis_size; print_idx++) {
            if (print_pretty_min_max_n_kstates.second <= print_idx) {
                std::cout << print_outer_prefix << message_prefix << data_tag
                          << "Break pretty print loop as the requested max number of contributions are already printed."
                          << std::endl;
                break;
            }
            unsigned idx = kstateIdx_and_probability_vector[print_idx].idx;
            const auto kstate_ptr = basis.vec_index()[idx];
            const std::complex<double> amplitude = eigien_vector(idx);
            const double probability = std::norm(amplitude);
            accumulated_probability += probability;
            if (print_pretty_min_max_n_kstates.first <= print_idx &&  probability < print_pretty_probability_treshold) {
                std::cout << print_outer_prefix << message_prefix << data_tag
                          <<  "Break pretty print loop as the rest of the kstate contributions have probabilities lower than the requested threshold." << std::endl;
                break;
            }
            std::cout << print_outer_prefix << message_prefix << data_tag
                      << "eigen vector no: " << std::setw(6) << n_eigien_vector << ", "
                      << "kstate constribution: " << (*kstate_ptr) << ", "
                      << "probability: " << std::noshowpos << probability << ", "
                      << "amplitude: " << std::showpos << presentation_layer_amplitude_factor * amplitude << ", "
                      << "kstate basis idx: " << std::setw(6) << idx << ", "
                      << "kstate print idx: " << std::setw(6) << print_idx
                      << "." << std::endl;
        } // end of print_idx loop
        std::cout << print_outer_prefix << message_prefix << data_tag
                  << "The listed above kstate contributions accumulated_probability: " << std::noshowpos << accumulated_probability
                  << "." << std::endl;
    } // end of n_eigien_vector loop
}

}

// #######################################################################
// ## do_common_recipe                                                  ##
// #######################################################################

namespace bfpt_common {

struct CommonRecipeReceipt {
    double energy;
    std::optional<arma::cx_vec> eigen_vector;
};

template<typename KstateT>
utility::Result<CommonRecipeReceipt, std::runtime_error>
do_common_recipe(const IKstateBasisPopulator<KstateT>& bais_populator,
                 const IKstateOperatorMatrix<KstateT>& hamiltonian,
                 kstate::Basis<KstateT>& basis,
                 const unsigned max_pt_order,
                 const unsigned k_n,
                 CommonRecipePrintFlags print_flags,
                 std::string print_outer_prefix = "",
                 unsigned n_threads = 1) {
    using ResultT = utility::Result<CommonRecipeReceipt, std::runtime_error>;
    using SiteStateT = typename KstateT::SiteType;
    assert(n_threads != 0);
    assert(n_threads <= 256);
    [[maybe_unused]] const size_t n_sites = basis.n_sites();
    assert(k_n < n_sites);
    arma::wall_clock timer;
    // --------------------------------------------------
    if (print_flags.print_unpopulated_basis_flag) {
        std::cout << print_outer_prefix << message_prefix << data_tag << "Unpopulated basis (0'th pt-order basis):" << std::endl;
        std::cout << basis;
    }
    if (print_flags.print_unpopulated_basis_size_flag) {
        std::cout << print_outer_prefix << message_prefix << data_tag << "Unpopulated basis (0'th pt-order basis) size: "
                  << basis.size() << "." << std::endl;
    }
    // --------------------------------------------------
    // Define hamiltonian and basis:
    std::cout << print_outer_prefix << message_prefix << progress_tag << "About to populate pt-basis." << std::endl;
    // Generate higher pt-orders subspace basis:
    timer.tic();
    generate_populated_basis(bais_populator, max_pt_order, basis, n_threads);
    const double time_populating_pt_basis = timer.toc();
    std::cout << print_outer_prefix << message_prefix << time_tag << "Populating pt-basis took " << time_populating_pt_basis << "s." << std::endl;
    std::cout << print_outer_prefix << message_prefix << progress_tag << "Has populated pt-basis." << std::endl;
    // --------------------------------------------------
    if (print_flags.print_populated_basis_flag) {
        std::cout << print_outer_prefix << message_prefix << data_tag << "Populated basis (" << max_pt_order << "'th pt-order basis):" << std::endl;
        std::cout << basis;
    }
    if (print_flags.print_populated_basis_size_flag) {
        std::cout << print_outer_prefix << message_prefix << data_tag << "Populated basis (" << max_pt_order << "'th pt-order basis) size: "
                  << basis.size() << "." << std::endl;
    }
    // --------------------------------------------------
    // Generate hamiltoniam matrix:
    std::cout << print_outer_prefix << message_prefix << progress_tag << "About to generate hamiltoniam." << std::endl;
    timer.tic();
    const auto kn_hamiltonian_matrix = generate_operator_matrix(hamiltonian, basis, k_n, n_threads);
    const double time_generating_kn_hamiltonian_matrix = timer.toc();
    std::cout << print_outer_prefix << message_prefix << time_tag << "Generating kn-hamiltoniam matrix took " << time_generating_kn_hamiltonian_matrix << "s." << std::endl;
    std::cout << print_outer_prefix << message_prefix << progress_tag << "Has generated kn-hamiltoniam." << std::endl;
    // --------------------------------------------------
    if (print_flags.print_hamiltonian_stats) {
        std::cout << print_outer_prefix << message_prefix << data_tag << "kn_hamiltonian_matrix.n_rows: " << kn_hamiltonian_matrix.n_rows << std::endl;
        std::cout << print_outer_prefix << message_prefix << data_tag << "kn_hamiltonian_matrix.n_cols: " << kn_hamiltonian_matrix.n_cols << std::endl;
        std::cout << print_outer_prefix << message_prefix << data_tag << "kn_hamiltonian_matrix.n_nonzero: " << kn_hamiltonian_matrix.n_nonzero << std::endl;
        const double fill_quotient = double(kn_hamiltonian_matrix.n_nonzero) / double(kn_hamiltonian_matrix.n_cols * kn_hamiltonian_matrix.n_rows);
        std::cout << print_outer_prefix << message_prefix << data_tag << "kn_hamiltonian_matrix.fill_quotient: " << fill_quotient * 100 << "%" << std::endl;
    }
    if (print_flags.print_sp_hamiltonian_flag) {
        std::cout << print_outer_prefix << message_prefix << data_tag << "kn_hamiltonian_matrix:";
        std::cout << kn_hamiltonian_matrix;
    }
    if (print_flags.print_hamiltonian_flag) {
        std::cout << print_outer_prefix << message_prefix << data_tag << "kn_hamiltonian_matrix:";
        std::cout << arma::cx_mat(kn_hamiltonian_matrix);
    }
    // --------------------------------------------------
    const bool should_calculate_eigenvectors = print_flags.print_eigen_vectors_flag || print_flags.print_pretty_vectors_flag || print_flags.print_density_operator_flag;
    // --------------------------------------------------
    if (should_calculate_eigenvectors) {
        std::cout << print_outer_prefix << message_prefix << progress_tag << "About to solve eigen problem (eigenvalues & eigenvectors)." << std::endl;
        timer.tic();
        const auto& eigs_sym_result = lin_alg::fallbacked_eigs_sym(lin_alg::WithVectors{}, kn_hamiltonian_matrix, 1, 1e-6);
        const double time_solving_eigen_problem = timer.toc();
        std::cout << print_outer_prefix << message_prefix << time_tag << "Solving eigen problem took: " << time_solving_eigen_problem << "s." << std::endl;
        if (eigs_sym_result.is_err()) {
            std::cout << print_outer_prefix << message_prefix << progress_tag << "Failed to solve eigen problem (eigenvalues & eigenvectors)." << std::endl;
            std::cout << print_outer_prefix << message_prefix << progress_tag << "The reported error message:" << std::endl;
            std::cout << eigs_sym_result.unwrap_err().what() << std::endl;
            return ResultT::Err(std::runtime_error(eigs_sym_result.unwrap_err().what()));
        }
        std::cout << print_outer_prefix << message_prefix << progress_tag << "Has solved eigen problem (eigenvalues & eigenvectors)." << std::endl;
        const auto eigen_info = eigs_sym_result.unwrap();
        const arma::vec& eigen_values = eigen_info.eigen_values;
        const arma::cx_mat& eigen_vectors = eigen_info.eigen_vectors;
        // --------------------------------------------------
        if (print_flags.print_eigen_values_flag) {
            std::cout << print_outer_prefix << message_prefix << data_tag << "eigen_values:" << std::endl;
            std::cout << print_outer_prefix << message_prefix << eigen_values;
        }
        if (print_flags.print_eigen_vectors_flag) {
            std::cout << print_outer_prefix << message_prefix << data_tag << "eigen_vectors:" << std::endl;
            std::cout << print_outer_prefix << message_prefix << eigen_vectors;
        }
        if (print_flags.print_pretty_vectors_flag) {
            pretty_print(basis, eigen_values, eigen_vectors,
                         print_flags.print_pretty_min_max_n_kstates, print_flags.print_pretty_probability_treshold,
                         print_outer_prefix);
        }
        // --------------------------------------------------
        if (print_flags.print_density_operator_flag) {
            const DensityOperator12<SiteStateT> density_operator_12 = calculate_reduced_density_operator_12<KstateT>(basis, eigen_vectors.col(0));
            pretty_density_operator_12(density_operator_12, print_outer_prefix);
        }
        // --------------------------------------------------
        return ResultT::Ok({eigen_values(0), eigen_vectors.col(0)});
    } else {
        std::cout << print_outer_prefix << message_prefix << progress_tag << "About to solve eigen problem (eigenvalues only)." << std::endl;
        timer.tic();
        const auto& eigs_sym_result = lin_alg::fallbacked_eigs_sym(lin_alg::WithoutVectors{}, kn_hamiltonian_matrix, 1, 1e-6);
        const double time_solving_eigen_problem = timer.toc();
        std::cout << print_outer_prefix << message_prefix << time_tag << "Solving eigen problem took: " << time_solving_eigen_problem << "s." << std::endl;
        if (eigs_sym_result.is_err()) {
            std::cout << print_outer_prefix << message_prefix << progress_tag << "Failed to solve eigen problem (eigenvalues only)." << std::endl;
            std::cout << print_outer_prefix << message_prefix << progress_tag << "The reported error message:" << std::endl;
            std::cout << eigs_sym_result.unwrap_err().what() << std::endl;
            return ResultT::Err(std::runtime_error(eigs_sym_result.unwrap_err().what()));
        }
        std::cout << print_outer_prefix << message_prefix << progress_tag << "Has solved eigen problem (eigenvalues only)." << std::endl;
        const arma::vec& eigen_values = eigs_sym_result.unwrap();
        // --------------------------------------------------
        if (print_flags.print_eigen_values_flag) {
            std::cout << print_outer_prefix << message_prefix << data_tag << "eigen_values:" << std::endl;
            std::cout << print_outer_prefix << message_prefix << eigen_values;
        }
        return ResultT::Ok({eigen_values(0), std::nullopt});
    }
    // --------------------------------------------------
    assert(false);
}

}  // namespace bfpt_common

#endif
