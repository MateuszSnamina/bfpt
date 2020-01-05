#include <model_monostar/lin_alg.hpp>

#include <cassert>

// #######################################################################
// ## re_to_cx and re_to_cx                                             ##
// #######################################################################

namespace lin_alg {

// arma::vec cx_to_re(const arma::cx_vec& v_cx) {
//     const auto size = v_cx.n_rows;
//     assert(size > 0);
//     arma::vec v_re(size * 2);
//     const arma::span span_re(0, size - 1);
//     const arma::span span_im(size, 2 * size - 1);
//     v_re(span_re) = arma::real(v_cx);
//     v_re(span_im) = arma::imag(v_cx);
//     return v_re;
// }

arma::cx_vec re_to_cx(const arma::vec& v_re) {
    assert(v_re.n_rows % 2 == 0);
    const auto size = v_re.n_rows / 2;
    assert(size > 0);
    arma::cx_vec v_cx(size);
    const arma::span span_re(0, size - 1);
    const arma::span span_im(size, 2 * size - 1);
    v_cx.set_real(v_re(span_re));
    v_cx.set_imag(v_re(span_im));
    return v_cx;
}

arma::cx_mat re_to_cx(const arma::mat& m_re) {
    assert(m_re.n_rows % 2 == 0);
    const auto size = m_re.n_rows / 2;
    assert(size > 0);
    arma::cx_mat m_cx(size, m_re.n_cols);
    const arma::span span_re(0, size - 1);
    const arma::span span_im(size, 2 * size - 1);
    m_cx.set_real(m_re.rows(span_re));
    m_cx.set_imag(m_re.rows(span_im));
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

}  // namespace lin_alg

// #######################################################################
// ## reduce_eigen_values                                               ##
// #######################################################################

namespace lin_alg {

arma::vec reduce_eigen_values(const arma::vec& eigen_values_not_reduced, double eps) {
    assert(eigen_values_not_reduced.n_rows % 2 == 0);
    arma::vec eigen_values = arma::vec(eigen_values_not_reduced.n_rows / 2);
    for (arma::uword i = 0; i < eigen_values_not_reduced.n_rows; i += 2) {
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
    for (const auto& span : spans) {
        std::cout << "     Degenracy subspace: "
                  << "[" << span.first << ", " << span.second << ")" << std::endl;
    }
    return spans;
}

}  // namespace lin_alg

// #######################################################################
// ## eig_sym                                                           ##
// #######################################################################

namespace lin_alg {

bool eig_sym(arma::vec& eigen_values, arma::cx_mat& eigen_vectors, const arma::cx_mat& matrix) {
    const arma::mat re_matrix = main_matrix_cx_to_re(matrix);
    arma::vec eigen_values_not_reduced;
    arma::mat re_eigen_vectors_not_reduced;
    const bool res = arma::eig_sym(eigen_values_not_reduced, re_eigen_vectors_not_reduced, re_matrix);
    const arma::cx_mat cx_eigen_vectors_not_reduced = re_to_cx(re_eigen_vectors_not_reduced);
    // Poor threshold hardcoding:
    const double eps = 1e-6;
    // Reduction --  eigen_values
    eigen_values = reduce_eigen_values(eigen_values_not_reduced, eps);
    // Analyse degeneracy subspaces:
    const std::vector<MySpan> spans = make_degeneracy_subspaces_analyse(eigen_values, eps);
    // Reduction --  eigen_vectors
    //eigen_vectors = cx_eigen_vectors_not_reduced;
    eigen_vectors = arma::cx_mat(matrix.n_cols, eigen_values.n_rows);
    for (unsigned span_idx = 0; span_idx < spans.size(); span_idx++) {
        const auto& span = spans[span_idx];
        assert(span.first < span.second);  // span cannot be empty.
        const arma::uword span_size = span.second - span.first;
        const arma::span not_reduced_span(2 * span.first, 2 * span.second - 1);
        const arma::span reduced_span(span.first, span.second - 1);
        arma::cx_mat basis = arma::orth(cx_eigen_vectors_not_reduced.cols(not_reduced_span));
        if (basis.n_cols != span_size) {
            std::cerr << "[warning] Warning for degeneracy subspace no " << span_idx << "(out of " << spans.size() << ")." << std::endl;
            std::cerr << "[warning] number of eigen_vectors was not reduced by factor of two." << std::endl;
            std::cerr << "[warning] this indicates the error, expect the subsapce is the last subspace" << std::endl;
            std::cerr << "[warning] In the latter case it means the subspace is not determined completely" << std::endl;
        }
        //TODO handle the error.
        eigen_vectors.cols(reduced_span) = basis;
    }
    return res;
}

}  // namespace lin_alg