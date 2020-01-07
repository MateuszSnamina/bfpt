#include <linear_algebra/linear_algebra.hpp>

#include <cassert>

// #######################################################################
// ## re_to_cx and re_to_cx                                             ##
// #######################################################################

namespace lin_alg {

// arma::vec cx_to_re(const arma::cx_vec& v_cx) {
//     assert(v_cx.n_rows > 0);
//     const auto size = v_cx.n_rows;
//     arma::vec v_re(size * 2);
//     const arma::span span_re(0, size - 1);
//     const arma::span span_im(size, 2 * size - 1);
//     v_re(span_re) = arma::real(v_cx);
//     v_re(span_im) = arma::imag(v_cx);
//     return v_re;
// }

arma::cx_vec re_to_cx(const arma::vec& v_re) {
    assert(v_re.n_rows % 2 == 0);
    assert(v_re.n_rows > 0);
    const auto size = v_re.n_rows / 2;
    arma::cx_vec v_cx(size);
    const arma::span span_re(0, size - 1);
    const arma::span span_im(size, 2 * size - 1);
    v_cx.set_real(v_re(span_re));
    v_cx.set_imag(v_re(span_im));
    return v_cx;
}

arma::sp_cx_vec re_to_cx(const arma::sp_vec& v_re) {
    assert(v_re.n_rows % 2 == 0);
    assert(v_re.n_rows > 0);
    const auto size = v_re.n_rows / 2;
    arma::cx_vec v_cx(size);
    using IterT = arma::sp_mat::const_iterator;
    for (IterT it = v_re.begin(); it != v_re.end(); ++it) {
        const auto& row = it.row();
        const auto& val = *it;
        const auto where_to_put = row < size ? row : row - size;
        const auto what_to_put = row < size ? std::complex<double>(val, 0) : std::complex<double>(0, val);
        v_cx(where_to_put) += what_to_put;
    }
    return v_cx;
}

arma::cx_mat re_to_cx(const arma::mat& m_re) {
    assert(m_re.n_rows % 2 == 0);
    assert(m_re.n_rows > 0);
    const auto size = m_re.n_rows / 2;
    arma::cx_mat m_cx(size, m_re.n_cols);
    const arma::span span_re(0, size - 1);
    const arma::span span_im(size, 2 * size - 1);
    m_cx.set_real(m_re.rows(span_re));
    m_cx.set_imag(m_re.rows(span_im));
    return m_cx;
}

arma::sp_cx_mat re_to_cx(const arma::sp_mat& m_re) {
    assert(m_re.n_rows % 2 == 0);
    assert(m_re.n_rows > 0);
    const auto size = m_re.n_rows / 2;
    arma::sp_cx_mat m_cx(size, m_re.n_cols);
    using IterT = arma::sp_mat::const_iterator;
    for (IterT it = m_re.begin(); it != m_re.end(); ++it) {
        const auto& col = it.col();
        const auto& row = it.row();
        const auto& val = *it;
        const auto where_to_put = row < size ? row : row - size;
        const auto what_to_put = row < size ? std::complex<double>(val, 0) : std::complex<double>(0, val);
        m_cx(where_to_put, col) += what_to_put;
    }
    return m_cx;
}

}  // namespace lin_alg

// #######################################################################
// ## main_matrix_cx_to_re                                              ##
// #######################################################################

namespace lin_alg {

arma::mat main_matrix_cx_to_re(const arma::cx_mat& m_cx) {
    assert(m_cx.n_rows == m_cx.n_cols);
    assert(m_cx.n_rows > 0);
    const auto size = m_cx.n_rows;
    arma::mat m_re(size * 2, size * 2);
    const arma::span span_re(0, size - 1);
    const arma::span span_im(size, 2 * size - 1);
    m_re(span_re, span_re) = +arma::real(m_cx);
    m_re(span_im, span_re) = +arma::imag(m_cx);
    m_re(span_re, span_im) = -arma::imag(m_cx);
    m_re(span_im, span_im) = +arma::real(m_cx);
    return m_re;
}

arma::sp_mat main_matrix_cx_to_re(const arma::sp_cx_mat& m_cx) {
    assert(m_cx.n_rows == m_cx.n_cols);
    assert(m_cx.n_rows > 0);
    const auto size = m_cx.n_rows;
    arma::sp_mat m_re(size * 2, size * 2);
    // const arma::span span_re(0, size - 1);
    // const arma::span span_im(size, 2 * size - 1);
    using IterT = arma::sp_cx_mat::const_iterator;
    for (IterT it = m_cx.begin(); it != m_cx.end(); ++it) {
        const auto& col = it.col();
        const auto& row = it.row();
        const auto& val = *it;
        const auto& val_real = val.real();
        const auto& val_imag = val.imag();
        m_re(row, col) = +val_real;
        m_re(row + size, col) = +val_imag;
        m_re(row, col + size) = -val_imag;
        m_re(row + size, col + size) = +val_real;
    }
    return m_re;
}

}  // namespace lin_alg

// #######################################################################
// ## reduce_eigen_values                                               ##
// #######################################################################

namespace lin_alg {

LinearAlgebraResult<arma::vec> reduce_eigen_values(const arma::vec& eigen_values_not_reduced, double eps) {
    assert(eigen_values_not_reduced.n_rows % 2 == 0);
    arma::vec eigen_values = arma::vec(eigen_values_not_reduced.n_rows / 2);
    for (arma::uword i = 0; i < eigen_values_not_reduced.n_rows; i += 2) {
        if (std::abs(eigen_values_not_reduced(i) - eigen_values_not_reduced(i + 1)) >= eps) {
            const std::string message1 = "The eigenvalues of the corresponding real eigen-problem do not follow the expected pattern.";
            const std::string message2 = "The following values are not 'near': eigen_values_not_reduced(i), eigen_values_not_reduced(i+1): " + std::to_string(eigen_values_not_reduced(i)) + ", " + std::to_string(eigen_values_not_reduced(i + 1)) + " (for i: " + std::to_string(i) + ").";
            const std::string message3 = "The 'near' relationship is consider with eps: " + std::to_string(eps) + ".";
            const std::string message = message1 + " " + message2 + " " + message3;
            std::cerr << "[debug-info] [line 1/3] " << message1 << std::endl;
            std::cerr << "[debug-info] [line 2/3] " << message2 << std::endl;
            std::cerr << "[debug-info] [line 3/3] " << message3 << std::endl;
            ReduceEigenValuesError details{i, eigen_values_not_reduced(i), eigen_values_not_reduced(i + 1), eps};
            return LinearAlgebraRuntimeException{message, details};
        }
        assert(std::abs(eigen_values_not_reduced(i) - eigen_values_not_reduced(i + 1)) < eps);
        eigen_values(i / 2) = eigen_values_not_reduced(i);
    }
    return eigen_values;
}

}  // namespace lin_alg

// #######################################################################
// ## make_degeneracy_subspaces_analyse                                 ##
// #######################################################################

namespace lin_alg {

std::vector<MySpan> make_degeneracy_subspaces_analyse(const arma::vec& eigen_values, double eps) {
    std::vector<MySpan> spans;
    arma::uword current_span_begin = 0;
    double current_value = eigen_values(0);
    for (arma::uword i = 1; i < eigen_values.n_rows; i++) {
        if (eigen_values(i) - current_value > eps) {
            spans.emplace_back(current_span_begin, i);
            current_span_begin = i;
            current_value = eigen_values(i);
        }
    }
    spans.emplace_back(current_span_begin, eigen_values.n_rows);
    // for (const auto& span : spans) {
    //     std::cout << "     Degenracy subspace: "
    //               << "[" << span.first << ", " << span.second << ")" << std::endl;
    // }
    return spans;
}

}  // namespace lin_alg

// #######################################################################
// ## eigs_sym                                                          ##
// #######################################################################

namespace lin_alg {

LinearAlgebraResult<HermitianEigenInfo>
eigs_sym(const arma::sp_cx_mat& matrix, unsigned n_vectors,
         unsigned n_extra_vectors, const char* form, double tol) {
    // --------------------------------------------------------------
    assert(n_vectors > 0);
    assert(matrix.n_cols == matrix.n_rows);
    const auto size = matrix.n_cols;
    const arma::sp_mat re_matrix = main_matrix_cx_to_re(matrix);
    assert(re_matrix.n_cols == re_matrix.n_rows);
    assert(re_matrix.n_cols == 2 * size);
    // --------------------------------------------------------------
    arma::vec eigen_values_not_reduced;
    arma::mat re_eigen_vectors_not_reduced;
    const arma::uword n_vectors_not_reduced_to_calculate = 2 * n_vectors + n_extra_vectors;
    assert(n_vectors_not_reduced_to_calculate < 2 * size);
    const bool res = arma::eigs_sym(eigen_values_not_reduced, re_eigen_vectors_not_reduced, re_matrix,
                                    n_vectors_not_reduced_to_calculate, form, tol);
    if (!res) {
        std::string message = "arma::eigs_sym claims it failed.";
        std::cerr << "[debug-info] " << message << "." << std::endl;
        return LinearAlgebraRuntimeException{message, ArmaEigsSymClaimsFailed{}};
    }
    if (eigen_values_not_reduced.n_rows < 2 * n_vectors) {
        const std::string message1 = "arma::eigs_sym does not claim it failed, but the number of eigenvalues it found is less than needed.";
        const std::string message2a = "(Requested " + std::to_string(n_vectors_not_reduced_to_calculate) + " eigenvalues, ";
        const std::string message2b = "got only " + std::to_string(eigen_values_not_reduced.n_rows) + " eigenvalues, ";
        const std::string message2c = "while needed " + std::to_string(2 * n_vectors) + " eigenvalues.)";
        const std::string message2 = message2a + message2b + message2c;
        const std::string message = message1 + message2;
        std::cerr << "[debug-info] [line 1/2] " << message1 << std::endl;
        std::cerr << "[debug-info] [line 2/2] " << message2 << std::endl;
        ArmaEigsSymFailedToFoundEnoughEigenvectors error_detail{n_vectors_not_reduced_to_calculate, eigen_values_not_reduced.n_rows, 2 * n_vectors};
        return LinearAlgebraRuntimeException{message, error_detail};
    }
    assert(re_eigen_vectors_not_reduced.n_rows == 2 * size);
    if (re_eigen_vectors_not_reduced.n_cols < 2 * n_vectors) {
        const std::string message1 = "arma::eigs_sym does not claim it failed, but the number of eigenvectors it found is less than needed.";
        const std::string message2a = "(Requested: " + std::to_string(n_vectors_not_reduced_to_calculate) + " eigenvectors, ";
        const std::string message2b = "got only: " + std::to_string(re_eigen_vectors_not_reduced.n_rows) + " eigenvectors, ";
        const std::string message2c = "while needed: " + std::to_string(2 * n_vectors) + " eigenvectors.)";
        const std::string message2 = message2a + message2b + message2c;
        const std::string message = message1 + " " + message2;
        std::cerr << "[debug-info] [line 1/2] " << message1 << std::endl;
        std::cerr << "[debug-info] [line 2/2] " << message2 << std::endl;
        ArmaEigsSymFailedToFoundEnoughEigenvectors error_detail{n_vectors_not_reduced_to_calculate, eigen_values_not_reduced.n_rows, 2 * n_vectors};
        return LinearAlgebraRuntimeException{message, error_detail};
        // return false;
    }
    eigen_values_not_reduced = eigen_values_not_reduced(arma::span(0, 2 * n_vectors - 1));
    re_eigen_vectors_not_reduced = re_eigen_vectors_not_reduced.cols(arma::span(0, 2 * n_vectors - 1));
    const arma::cx_mat cx_eigen_vectors_not_reduced = re_to_cx(re_eigen_vectors_not_reduced);
    // --------------------------------------------------------------
    // Reduction - eigen_values:
    const auto reduce_eigen_values_result = reduce_eigen_values(eigen_values_not_reduced, 1000 * tol);
    if (reduce_eigen_values_result.is_err()) {
        return reduce_eigen_values_result.unwrap_err();
    }
    const arma::vec eigen_values = reduce_eigen_values_result.unwrap();
    // --------------------------------------------------------------
    // Reduction - eigen_vectors:
    const std::vector<MySpan> spans = make_degeneracy_subspaces_analyse(eigen_values, tol);
    arma::cx_mat eigen_vectors(matrix.n_cols, eigen_values.n_rows);
    for (unsigned span_idx = 0; span_idx < spans.size(); span_idx++) {
        const auto& span = spans[span_idx];
        assert(span.first < span.second);  // span cannot be empty.
        const arma::uword span_size = span.second - span.first;
        const arma::span not_reduced_span(2 * span.first, 2 * span.second - 1);
        const arma::span reduced_span(span.first, span.second - 1);
        arma::cx_mat basis = arma::orth(cx_eigen_vectors_not_reduced.cols(not_reduced_span), 100 * tol);
        const bool got_less_than_expected = basis.n_cols < span_size;
        const bool got_more_than_expected = basis.n_cols > span_size;
        const bool is_the_last_subspace = span_idx + 1 == spans.size();
        if (got_less_than_expected || (got_more_than_expected && !is_the_last_subspace)) {
            const std::string message1 = "Failed to reproduce the basis of a complex degeneracy subspace.";
            const std::string message2a = "(Expected the complex degeneracy subspace dimension: " + std::to_string(span_size) + ", ";
            const std::string message2b = "got complex subspace degeneracy of dimension: " + std::to_string(basis.n_cols) + ".)";
            const std::string message3 = "The error occured for degeneracy subspace no " + std::to_string(span_idx) + " (out of " + std::to_string(spans.size()) + ").";
            const std::string message2 = message2a + message2b;
            const std::string message = message1 + " " + message2 + " " + message3;
            std::cerr << "[debug-info] [line 1/4] " << message1 << std::endl;
            std::cerr << "[debug-info] [line 2/4] " << message2a << std::endl;
            std::cerr << "[debug-info] [line 3/4] " << message2b << std::endl;
            std::cerr << "[debug-info] [line 4/4] " << message3 << std::endl;
            FailedToReproduceComplexDegeneracySubspace error_detail{span, span_size, basis.n_cols};
            return LinearAlgebraRuntimeException{message, error_detail};
        }
        eigen_vectors.cols(reduced_span) = basis;
    }
    return HermitianEigenInfo{eigen_values, eigen_vectors};
}

}  // namespace lin_alg