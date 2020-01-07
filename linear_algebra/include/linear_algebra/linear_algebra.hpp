
#ifndef LINEAR_ALGEBRA_LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_LINEAR_ALGEBRA_HPP

#include <linear_algebra/result.hpp>

#include <armadillo>

#include <stdexcept>
#include <utility>
#include <vector>

// #######################################################################
// ## Types:                                                            ##
// #######################################################################

namespace lin_alg {

/*
 * The type for inclusive-exclusive span definition.
 * (as opossed to arma::span being inclusive-inclusive).
 */
using MySpan = std::pair<arma::uword, arma::uword>;

/*
 * Typical return type wrapper.
 */
template <typename OutT>
using LinearAlgebraResult = Result<OutT, LinearAlgebraRuntimeException>;

// #######################################################################
// ## Functions:                                                        ##
// #######################################################################

/*
 * Functions mapping cx_mat (or cx_vec)
 * into mat (or vec) with two times more elements.
 * The mapping is performed column-wise, so that
 * each column in real versions has the following block layout:
 * ╭                        ╮
 * │ Re(corrspondin_cx_col) │
 * │ Im(corrspondin_cx_col) │
 * ╰                        ╯
 * 
 * Example (of a cx_vec-vec mapping):
 * cx_vec{1+2i,3+4i}  <-mapping->  vec{1,3,2,4}
 * 
 * Example (of a cx_vec-vec mapping):
 * ╭      ╮               ╭   ╮
 * │ 1+2i │  <-mapping->  │ 1 │
 * │ 3+4i │               │ 3 │
 * ╰      ╯               │ 2 │
 *                        │ 4 │
 *                        ╰   ╯
 * Example (of a cx_mat-mat mapping):
 * ╭              ╮               ╭       ╮
 * │ 1+2i, 11+12i │  <-mapping->  │ 1, 11 │
 * │ 3+4i, 13+14i │               │ 3, 13 │
 * ╰              ╯               │ 2, 12 │
 *                                │ 4, 14 │
 *                                ╰       ╯
 */
// arma::vec cx_to_re(const arma::cx_vec& v_cx);
arma::cx_vec re_to_cx(const arma::vec& v_re);
arma::sp_cx_vec re_to_cx(const arma::sp_vec& v_re);
arma::cx_mat re_to_cx(const arma::mat& m_re);
arma::sp_cx_mat re_to_cx(const arma::sp_mat& m_re);

/*
 * Function making mat from cx_mat,
 * so that the result matrix block layout is as follow:
 * ╭                       ╮
 * │  +Re(m_cx), -Im(m_cx) │
 * │  +Im(m_cx), +Re(m_cx) │
 * ╰                       ╯
 * where `m_cx` is the imput matrix of `cx_mat` type.
 * 
 * Note: If the input is a hermitian matrix,
 *       then the result matrix is symmetric.
 */
arma::mat main_matrix_cx_to_re(const arma::cx_mat& m_cx);

/*
 * Function eliminating every second element in a vec.
 * ╭   ╮                   ╭   ╮
 * │ 1 │  --function-->    │ 1 │
 * │ 1 │                   │ 3 │
 * │ 3 │                   ╰   ╯
 * │ 3 │
 * ╰   ╯
 * 
 * The input vector is assumed to be filled in a way that:
 * the n-th and the (n+1)-th have the save value
 * (do not differ more than eps).
 * The function asserts the confition is fulfilled.
 */

struct ReduceEigenValuesError {
    arma::uword i;        // index for eigen_values_not_reduced
    double value_at_i;    // value at eigen_values_not_reduced(i)
    double value_at_ip1;  // value at eigen_values_not_reduced(i+1)
    double assumed_threshold;
};

LinearAlgebraResult<arma::vec> reduce_eigen_values(const arma::vec& eigen_values_not_reduced, double eps);

/*
 * Function making description of degeneracy_subspaces.
 * ╭   ╮                   ╭              ╮
 * │ 0 │  --function-->    │ MySpan(0, 1) │
 * │ 1 │                   │ MySpan(1, 4) │
 * │ 1 │                   │ MySpan(4, 6) │
 * │ 1 │                   ╰              ╯
 * │ 3 │
 * │ 3 │
 * ╰   ╯
 * 
 * The input eigen_values is assumed to be sorted (ascending).
 * The two eigen_vectors are treated as from the same degeneracy subspaces
 * if their eigen values differ no more than eps.
 */
std::vector<MySpan> make_degeneracy_subspaces_analyse(const arma::vec& eigen_values, double eps);

/*
 * My version of arma::eig_sym(...) function.
 * 
 * My version solves eingen problem for cx_mat
 * by translating it into eingen problem for mat. 
 */

bool eig_sym(arma::vec& eigen_values, const arma::cx_mat& matrix);
bool eig_sym(arma::vec& eigen_values, arma::cx_mat& eigen_vectors, const arma::cx_mat& matrix);
bool eigs_sym(arma::vec& eigen_values, const arma::sp_cx_mat& matrix,
              unsigned n_vectors, unsigned n_extra_vectors, const char* form, double tol);
bool eigs_sym(arma::vec& eigen_values, arma::cx_mat& eigen_vectors, const arma::sp_cx_mat& matrix,
              unsigned n_vectors, unsigned n_extra_vectors, const char* form, double tol);

struct HermitianEigenInfo {
    arma::vec eigen_values;
    arma::cx_mat eigen_vectors;
};

struct ArmaEigsSymClaimsFailed {
};

struct ArmaEigsSymFailedToFoundEnoughEigenvectors {
    arma::uword n_expected_number_of_eigenvectors;
    arma::uword n_got_number_of_eigenvectors;
    arma::uword n_needed_eigenvectors;
};

struct FailedToReproduceComplexDegeneracySubspace {
    MySpan span;
    arma::uword n_expected_degeneracy_subspace_dimension;
    arma::uword n_got_degeneracy_subspace_dimension;

};

LinearAlgebraResult<HermitianEigenInfo>
eigs_sym(const arma::sp_cx_mat& matrix, unsigned n_vectors,
              unsigned n_extra_vectors, const char* form, double tol);

}  // namespace lin_alg

#endif