#pragma once

#include <string>

namespace so_app {

struct RawProgramOptions {
    unsigned n_sites;
    unsigned n_pt;
    std::string n_max_site_spin_excitations_string;
    std::string n_max_site_orbit_excitations_string;
    bool accept_orbit_site_excitations_only_if_near_domain_wall;
    std::string model_type_string;
    double hamiltonian_s_coef;
    double hamiltonian_ss_coef;
    double hamiltonian_tau_z_coef;
    double hamiltonian_tau_minus_coef;
    double hamiltonian_Pzz_coef;
    double hamiltonian_Pxz_coef;
    double hamiltonian_Pxx_coef;
    double hamiltonian_ss_Pzz_coef;
    double hamiltonian_ss_Pxz_coef;
    double hamiltonian_ss_Pxx_coef;
    double hamiltonian_free_coef;
    std::string orbital_theta_string;
    double average_ss;  // value used for `auto`-type calculation of thera_opt
    std::string run_type_string;
    bool print_unpopulated_basis_flag;
    bool print_unpopulated_basis_size_flag;
    bool print_populated_basis_flag;
    bool print_populated_basis_size_flag;
    bool print_hamiltonian_stats;
    bool print_sp_hamiltonian_flag;
    bool print_hamiltonian_flag;
    bool print_eigen_values_flag;
    bool print_eigen_vectors_flag;
    bool print_pretty_vectors_flag;
    bool print_density_operator_flag;
    unsigned print_pretty_min_n_kstates;
    unsigned print_pretty_max_n_kstates;
    double print_pretty_probability_treshold;
    std::string es_momentum_domain_string;
    unsigned es_n_k;
    unsigned n_threads;
};

RawProgramOptions grep_program_options(int argc, char** argv);

}  // end of namespace so_app
