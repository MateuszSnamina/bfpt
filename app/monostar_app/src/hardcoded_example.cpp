
#include <monostar_app/hardcoded_example.hpp>

#include <armadillo>

#include <array>
#include <complex>

namespace monostar_app {

void do_hardcoded_example_analyse() {
    std::array<int, 8> v00 = {1, 0, 0, 0, 0, 0, 0, 0};
    std::array<int, 8> v01 = {0, 1, 0, 0, 0, 0, 0, 0};
    std::array<int, 8> v02 = {0, 0, 1, 0, 0, 0, 0, 0};
    std::array<int, 8> v03 = {0, 0, 0, 1, 0, 0, 0, 0};
    std::array<int, 8> v04 = {0, 0, 0, 0, 1, 0, 0, 0};
    std::array<int, 8> v05 = {0, 0, 0, 0, 0, 1, 0, 0};
    std::array<int, 8> v06 = {0, 0, 0, 0, 0, 0, 1, 0};
    std::array<int, 8> v07 = {0, 0, 0, 0, 0, 0, 0, 1};

    std::array<int, 8> v08 = {1, 1, 1, 0, 0, 0, 0, 0};
    std::array<int, 8> v09 = {0, 1, 1, 1, 0, 0, 0, 0};
    std::array<int, 8> v10 = {0, 0, 1, 1, 1, 0, 0, 0};
    std::array<int, 8> v11 = {0, 0, 0, 1, 1, 1, 0, 0};
    std::array<int, 8> v12 = {0, 0, 0, 0, 1, 1, 1, 0};
    std::array<int, 8> v13 = {0, 0, 0, 0, 0, 1, 1, 1};
    std::array<int, 8> v14 = {1, 0, 0, 0, 0, 0, 1, 1};
    std::array<int, 8> v15 = {1, 1, 0, 0, 0, 0, 0, 1};

    std::array<int, 8> v16 = {1, 1, 0, 0, 0, 0, 1, 0};
    std::array<int, 8> v17 = {0, 1, 1, 0, 0, 0, 0, 1};
    std::array<int, 8> v18 = {1, 0, 1, 1, 0, 0, 0, 0};
    std::array<int, 8> v19 = {0, 1, 0, 1, 1, 0, 0, 0};
    std::array<int, 8> v20 = {0, 0, 1, 0, 1, 1, 0, 0};
    std::array<int, 8> v21 = {0, 0, 0, 1, 0, 1, 1, 0};
    std::array<int, 8> v22 = {0, 0, 0, 0, 1, 0, 1, 1};
    std::array<int, 8> v23 = {1, 0, 0, 0, 0, 1, 0, 1};

    std::array<int, 8> v24 = {1, 1, 0, 0, 0, 1, 0, 0};
    std::array<int, 8> v25 = {0, 1, 1, 0, 0, 0, 1, 0};
    std::array<int, 8> v26 = {0, 0, 1, 1, 0, 0, 0, 1};
    std::array<int, 8> v27 = {1, 0, 0, 1, 1, 0, 0, 0};
    std::array<int, 8> v28 = {0, 1, 0, 0, 1, 1, 0, 0};
    std::array<int, 8> v29 = {0, 0, 1, 0, 0, 1, 1, 0};
    std::array<int, 8> v30 = {0, 0, 0, 1, 0, 0, 1, 1};
    std::array<int, 8> v31 = {1, 0, 0, 0, 1, 0, 0, 1};

    std::array<int, 8> v32 = {1, 1, 0, 0, 1, 0, 0, 0};
    std::array<int, 8> v33 = {0, 1, 1, 0, 0, 1, 0, 0};
    std::array<int, 8> v34 = {0, 0, 1, 1, 0, 0, 1, 0};
    std::array<int, 8> v35 = {0, 0, 0, 1, 1, 0, 0, 1};
    std::array<int, 8> v36 = {1, 0, 0, 0, 1, 1, 0, 0};
    std::array<int, 8> v37 = {0, 1, 0, 0, 0, 1, 1, 0};
    std::array<int, 8> v38 = {0, 0, 1, 0, 0, 0, 1, 1};
    std::array<int, 8> v39 = {1, 0, 0, 1, 0, 0, 0, 1};

    std::array<int, 8> v40 = {1, 1, 0, 1, 0, 0, 0, 0};
    std::array<int, 8> v41 = {0, 1, 1, 0, 1, 0, 0, 0};
    std::array<int, 8> v42 = {0, 0, 1, 1, 0, 1, 0, 0};
    std::array<int, 8> v43 = {0, 0, 0, 1, 1, 0, 1, 0};
    std::array<int, 8> v44 = {0, 0, 0, 0, 1, 1, 0, 1};
    std::array<int, 8> v45 = {1, 0, 0, 0, 0, 1, 1, 0};
    std::array<int, 8> v46 = {0, 1, 0, 0, 0, 0, 1, 1};
    std::array<int, 8> v47 = {1, 0, 1, 0, 0, 0, 0, 1};

    std::array<std::array<int, 8>, 48> basis = {
        v00,
        v01,
        v02,
        v03,
        v04,
        v05,
        v06,
        v07,
        v08,
        v09,
        v10,
        v11,
        v12,
        v13,
        v14,
        v15,
        v16,
        v17,
        v18,
        v19,
        v20,
        v21,
        v22,
        v23,
        v24,
        v25,
        v26,
        v27,
        v28,
        v29,
        v30,
        v31,
        v32,
        v33,
        v34,
        v35,
        v36,
        v37,
        v38,
        v39,
        v40,
        v41,
        v42,
        v43,
        v44,
        v45,
        v46,
        v47};
    arma::cx_mat H(48, 48, arma::fill::zeros);
    for (int i = 0; i < 16; i++) {
        H(i, i) = -1.;
    }
    for (int i = 0; i < 48; i++) {
        const std::array<int, 8> bra = basis[i];
        for (int j = 0; j < 48; j++) {  //DEBUG
            const std::array<int, 8> ket = basis[j];
            for (int h = 0; h < 8; h++) {
                const int hp1 = (h + 1) % 8;
                if (ket[h] == 0 && ket[hp1] == 0) {
                    std::array<int, 8> h_ket = ket;
                    h_ket[h] = 1;
                    h_ket[hp1] = 1;
                    if (h_ket == bra) {
                        H(i, j) = 0.5;
                    }
                }
                if (ket[h] == 1 && ket[hp1] == 1) {
                    std::array<int, 8> h_ket = ket;
                    h_ket[h] = 0;
                    h_ket[hp1] = 0;
                    if (h_ket == bra) {
                        H(i, j) = 0.5;
                    }
                }
            }
        }
    }
    // ----
    arma::vec eigen_values;
    arma::cx_mat eigen_vectors;
    arma::eig_sym(eigen_values, eigen_vectors, H);
    // ----
    std::cout << eigen_values;
    std::cout << "min eigen_value: " << eigen_values(0) << std::endl;
    std::cout << "eigvec2, poczatek*3: " << std::endl;
    std::cout << eigen_vectors(arma::span(0, 23), 2) << std::endl;
    std::cout << "eigvec2, poczatek: " << std::endl;
    std::cout << eigen_vectors(arma::span(0, 7), 2) << std::endl;
    std::cout << "eigvec3, poczatek: " << std::endl;
    std::cout << eigen_vectors(arma::span(0, 7), 3) << std::endl;
    std::cout << "eigvec4, poczatek: " << std::endl;
    std::cout << eigen_vectors(arma::span(0, 7), 4) << std::endl;
    std::cout << "eigvec5, poczatek: " << std::endl;
    std::cout << eigen_vectors(arma::span(0, 7), 5) << std::endl;
    // ----
    const double sq2 = std::sqrt(2);
    const std::complex<double> chi(1 / sq2, 1 / sq2);
    const std::complex<double> chi3(-1 / sq2, 1 / sq2);
    const std::complex<double> cxi(0, 1);
    arma::cx_vec8 t17 = {2. / 4, sq2 / 4, 0, -sq2 / 4, -2. / 4, -sq2 / 4, 0, sq2 / 4};
    arma::cx_vec8 t1735 = {1 / sq2, 0, 0, 0, -1 / sq2, 0, 0, 0};
    arma::cx_vec8 t0 = {1. / 2 / sq2, 1. / 2 / sq2, 1. / 2 / sq2, 1. / 2 / sq2, 1. / 2 / sq2, 1. / 2 / sq2, 1. / 2 / sq2, 1. / 2 / sq2};
    arma::cx_vec8 t4 = {1. / 2 / sq2, -1. / 2 / sq2, 1. / 2 / sq2, -1. / 2 / sq2, 1. / 2 / sq2, -1. / 2 / sq2, 1. / 2 / sq2, -1. / 2 / sq2};
    arma::cx_vec8 t1 = {1. / 2 / sq2, chi / 2.0 / sq2, cxi / 2.0 / sq2, chi3 / 2.0 / sq2, -1. / 2.0 / sq2, -chi / 2.0 / sq2, -cxi / 2.0 / sq2, -chi3 / 2.0 / sq2};
    arma::cx_vec8 t7 = {1. / 2 / sq2, -chi3 / 2.0 / sq2, -cxi / 2.0 / sq2, -chi / 2.0 / sq2, -1. / 2.0 / sq2, +chi3 / 2.0 / sq2, +cxi / 2.0 / sq2, +chi / 2.0 / sq2};
    arma::cx_vec8 t3 = {1. / 2 / sq2, +chi3 / 2.0 / sq2, -cxi / 2.0 / sq2, +chi / 2.0 / sq2, -1. / 2.0 / sq2, -chi3 / 2.0 / sq2, +cxi / 2.0 / sq2, -chi / 2.0 / sq2};
    // arma::vec8 t = {1/sq2, 0, 0, 0, -1/sq2, 0, 0, 0};
    arma::cx_vec8 t = t7;  //(t1 + t7) / sq2;
    arma::cx_mat beta(48, 6, arma::fill::zeros);
    beta(arma::span(0, 7), 0) = t;
    beta(arma::span(8, 15), 1) = t;
    beta(arma::span(16, 23), 2) = t;
    beta(arma::span(24, 31), 3) = t;
    beta(arma::span(32, 39), 4) = t;
    beta(arma::span(40, 47), 5) = t;
    std::cout << beta << std::endl;
    arma::cx_mat H1 = beta.t() * H * beta;
    // ----
    arma::vec eigen_values_1;
    arma::cx_mat eigen_vectors_1;
    arma::eig_sym(eigen_values_1, eigen_vectors_1, H1);
    // ----
    std::cout << "eigen_values_1" << std::endl;
    std::cout << eigen_values_1;
    std::cout << "eigen_vectors_1.col(0)" << std::endl;
    std::cout << eigen_vectors_1.col(0) / eigen_vectors_1(0, 0);
}

}  // namespace monostar_app
