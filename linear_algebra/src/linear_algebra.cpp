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
            const std::string message = "i, eigen_values_not_reduced(i), and eigen_values_not_reduced(i+1): " + std::to_string(i) + ", " + std::to_string(eigen_values_not_reduced(i)) + ", " + std::to_string(eigen_values_not_reduced(i + 1)) + ".";
            std::cerr << "[debug-info] " << message << std::endl;
            ReduceEigenValuesError details{i, eigen_values_not_reduced(i), eigen_values_not_reduced(i + 1), eps};
            LinearAlgebraRuntimeException error(message, details);
            return error;
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
// ## eig_sym                                                           ##
// #######################################################################

namespace lin_alg {

bool eig_sym(arma::vec& eigen_values, const arma::cx_mat& matrix) {
    // --------------------------------------------------------------
    assert(matrix.n_cols == matrix.n_rows);
    const double eps = 1e-6;
    const arma::mat re_matrix = main_matrix_cx_to_re(matrix);
    // --------------------------------------------------------------
    arma::vec eigen_values_not_reduced;
    const bool res = arma::eig_sym(eigen_values_not_reduced, re_matrix);
    if (!res) {
        std::cerr << "[warning] arma::eig_sym claims it failed." << std::endl;
        return false;
    }
    // --------------------------------------------------------------
    eigen_values = reduce_eigen_values(eigen_values_not_reduced, eps).unwrap();
    return true;
}

bool eig_sym(arma::vec& eigen_values, arma::cx_mat& eigen_vectors, const arma::cx_mat& matrix) {
    // --------------------------------------------------------------
    assert(matrix.n_cols == matrix.n_rows);
    const auto size = matrix.n_cols;
    const double eps = 1e-6;
    const arma::mat re_matrix = main_matrix_cx_to_re(matrix);
    // --------------------------------------------------------------
    arma::vec eigen_values_not_reduced;
    arma::mat re_eigen_vectors_not_reduced;
    const bool res = arma::eig_sym(eigen_values_not_reduced, re_eigen_vectors_not_reduced, re_matrix);
    if (!res) {
        std::cerr << "[warning] arma::eig_sym claims it failed." << std::endl;
        return false;
    }
    assert(eigen_values_not_reduced.n_rows == 2 * size);
    assert(re_eigen_vectors_not_reduced.n_rows == 2 * size);
    assert(re_eigen_vectors_not_reduced.n_cols == 2 * size);
    const arma::cx_mat cx_eigen_vectors_not_reduced = re_to_cx(re_eigen_vectors_not_reduced);
    // --------------------------------------------------------------
    // Reduction -- eigen_values:
    // //Future error handling:
    // const auto reduce_eigen_values_result = reduce_eigen_values(eigen_values_not_reduced, eps);
    // if (reduce_eigen_values_result.is_err()) {
    //   return TheFunctionResultType::Err(reduce_eigen_values_result.unwrap_err());
    // }
    // eigen_values = reduce_eigen_values_result.unwrap();
    eigen_values = reduce_eigen_values(eigen_values_not_reduced, eps).unwrap();

    // --------------------------------------------------------------
    // Analyse degeneracy subspaces:
    const std::vector<MySpan> spans = make_degeneracy_subspaces_analyse(eigen_values, eps);
    // Reduction - eigen_vectors:
    eigen_vectors = arma::cx_mat(size, eigen_values.n_rows);
    for (unsigned span_idx = 0; span_idx < spans.size(); span_idx++) {
        const auto& span = spans[span_idx];
        assert(span.first < span.second);  // span cannot be empty.
        const arma::uword span_size = span.second - span.first;
        const arma::span not_reduced_span(2 * span.first, 2 * span.second - 1);
        const arma::span reduced_span(span.first, span.second - 1);
        arma::cx_mat basis = arma::orth(cx_eigen_vectors_not_reduced.cols(not_reduced_span));
        if (basis.n_cols != span_size) {
            std::cerr << "[warning] Warning for degeneracy subspace no " << span_idx << " (out of " << spans.size() << ")." << std::endl;
            std::cerr << "[warning] Number of eigen-vectors was not reduced by factor of two." << std::endl;
            std::cerr << "[warning] Number of eigen-vectors before reduction: " << 2 * (span.second - span.first) << "." << std::endl;
            std::cerr << "[warning] Number of eigen-vectors after reduction: " << basis.n_cols << "." << std::endl;
            std::cerr << "[warning] This indicates an error." << std::endl;
            return false;
        }
        eigen_vectors.cols(reduced_span) = basis;
    }
    return true;
}

}  // namespace lin_alg

// #######################################################################
// ## eigs_sym                                                          ##
// #######################################################################

namespace lin_alg {

bool eigs_sym(arma::vec& eigen_values, const arma::sp_cx_mat& matrix,
              unsigned n_vectors, unsigned n_extra_vectors, const char* form, double tol) {
    // --------------------------------------------------------------
    assert(n_vectors > 0);
    assert(matrix.n_cols == matrix.n_rows);
    const auto size = matrix.n_cols;
    const arma::sp_mat re_matrix = main_matrix_cx_to_re(matrix);
    assert(re_matrix.n_cols == re_matrix.n_rows);
    assert(re_matrix.n_cols == 2 * size);
    // --------------------------------------------------------------
    arma::vec eigen_values_not_reduced;
    assert(2 * n_vectors + n_extra_vectors < 2 * size);
    const bool res = arma::eigs_sym(eigen_values_not_reduced, re_matrix,
                                    2 * n_vectors + n_extra_vectors, form, tol);
    if (!res) {
        std::cerr << "[warning] arma::eigs_sym claims it failed." << std::endl;
        return false;
    }
    if (eigen_values_not_reduced.n_rows < 2 * n_vectors) {
        std::cerr << "[warning] arma::eigs_sym does not claim it failed," << std::endl;
        std::cerr << "[warning] but the number of found eigen_values is too small." << std::endl;
        return false;
    }
    eigen_values_not_reduced = eigen_values_not_reduced(arma::span(0, 2 * n_vectors - 1));
    // --------------------------------------------------------------
    eigen_values = reduce_eigen_values(eigen_values_not_reduced, 100 * tol).unwrap();
    return true;
}

bool eigs_sym(arma::vec& eigen_values, arma::cx_mat& eigen_vectors, const arma::sp_cx_mat& matrix,
              unsigned n_vectors, unsigned n_extra_vectors, const char* form, double tol) {
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
    assert(2 * n_vectors + n_extra_vectors < 2 * size);
    const bool res = arma::eigs_sym(eigen_values_not_reduced, re_eigen_vectors_not_reduced, re_matrix,
                                    2 * n_vectors + n_extra_vectors, form, tol);
    if (!res) {
        std::cerr << "[warning] arma::eigs_sym claims it failed." << std::endl;
        return false;
    }
    if (eigen_values_not_reduced.n_rows < 2 * n_vectors) {
        std::cerr << "[warning] arma::eigs_sym does not claim it failed," << std::endl;
        std::cerr << "[warning] but the number of found eigen_values is too small." << std::endl;
        return false;
    }
    assert(re_eigen_vectors_not_reduced.n_rows == 2 * size);
    if (re_eigen_vectors_not_reduced.n_cols < 2 * n_vectors) {
        std::cerr << "[warning] arma::eigs_sym does not claim it failed," << std::endl;
        std::cerr << "[warning] but the number of found eigen_vectors is too small." << std::endl;
        return false;
    }
    eigen_values_not_reduced = eigen_values_not_reduced(arma::span(0, 2 * n_vectors - 1));
    re_eigen_vectors_not_reduced = re_eigen_vectors_not_reduced.cols(arma::span(0, 2 * n_vectors - 1));
    const arma::cx_mat cx_eigen_vectors_not_reduced = re_to_cx(re_eigen_vectors_not_reduced);
    // --------------------------------------------------------------
    // Reduction - eigen_values:
    eigen_values = reduce_eigen_values(eigen_values_not_reduced, 1000 * tol).unwrap();
    // --------------------------------------------------------------
    // Reduction - eigen_vectors:
    const std::vector<MySpan> spans = make_degeneracy_subspaces_analyse(eigen_values, tol);
    eigen_vectors = arma::cx_mat(matrix.n_cols, eigen_values.n_rows);
    for (unsigned span_idx = 0; span_idx < spans.size(); span_idx++) {
        const auto& span = spans[span_idx];
        assert(span.first < span.second);  // span cannot be empty.
        const arma::uword span_size = span.second - span.first;
        const arma::span not_reduced_span(2 * span.first, 2 * span.second - 1);
        const arma::span reduced_span(span.first, span.second - 1);
        arma::cx_mat basis = arma::orth(cx_eigen_vectors_not_reduced.cols(not_reduced_span), 100 * tol);
        if (basis.n_cols != span_size) {
            std::cerr << "[warning] Warning for degeneracy subspace no " << span_idx << " (out of " << spans.size() << ")." << std::endl;
            std::cerr << "[warning] Number of eigen-vectors was not reduced by factor of two." << std::endl;
            std::cerr << "[warning] Number of eigen-vectors before reduction: " << 2 * (span.second - span.first) << "." << std::endl;
            std::cerr << "[warning] Number of eigen-vectors after reduction: " << basis.n_cols << "." << std::endl;
            std::cerr << "[warning] This indicates an error, expect one specific case:" << std::endl;
            std::cerr << "[warning] the subsapce is the last subspace, and" << std::endl;
            std::cerr << "[warning] the reduction gives rise to more than (number of eigen-vectors before reduction)/2 vectors" << std::endl;
            std::cerr << "[warning] In the latter case it means the subspace is not determined completely" << std::endl;
            // TODO handle the error.
        }
        eigen_vectors.cols(reduced_span) = basis;
    }
    return true;
}

}  // namespace lin_alg