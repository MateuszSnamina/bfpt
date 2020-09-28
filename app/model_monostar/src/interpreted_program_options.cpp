#include <model_monostar/interpreted_program_options.hpp>

#include <model_monostar/raw_program_options.hpp>
#include <model_monostar/interpret_model_type_string.hpp>
#include <model_monostar/interpret_run_type_string.hpp>
#include <model_monostar/interpret_es_momentum_domain.hpp>

#include <stdexcept>

InterpretedProgramOptions interpret_program_options(const RawProgramOptions& raw_program_options) {
    InterpretedProgramOptions interpreted_program_options;
    interpreted_program_options.n_sites = raw_program_options.n_sites;
    interpreted_program_options.n_pt = raw_program_options.n_pt;
    if (const auto _ = interpret_model_type_string(raw_program_options.model_type_string)) {
        interpreted_program_options.model_type = _.unwrap();
    } else {
        const std::string message1 = "Problem with interpreting 'model_type_string' program option.";
        const std::string message = message1 + "\n" + _.unwrap_err().what();
        throw std::runtime_error(message);
    }
    interpreted_program_options.hamiltonian_af_fm_params = HamiltonianAfFmParams::Builder()
            .set_J_classical(raw_program_options.hamiltonian_J_classical)
            .set_J_quantum(raw_program_options.hamiltonian_J_quantum)
            .set_B(raw_program_options.hamiltonian_B)
            .build();
    interpreted_program_options.hamiltonian_fo_params = HamiltonianFoParams::Builder()
            .set_Pdelta_coef(raw_program_options.hamiltonian_Pdelta_coef)
            .set_Pzz_coef(raw_program_options.hamiltonian_Pzz_coef)
            .set_Pxz_coef(raw_program_options.hamiltonian_Pxz_coef)
            .set_Pxx_coef(raw_program_options.hamiltonian_Pxx_coef)
            .build();
    interpreted_program_options.theta_opt = raw_program_options.theta_opt;
    if (const auto _ = interpret_run_type_string(raw_program_options.run_type_string)) {
        interpreted_program_options.run_type = _.unwrap();
    } else {
        const std::string message1 = "Problem with interpreting 'model_type_string' program option.";
        const std::string message = message1 + "\n" + _.unwrap_err().what();
        throw std::runtime_error(message);
    }
    if (const auto _ = interpret_es_momentum_domain(raw_program_options.es_momentum_domain_string)) {
        const auto es_momentum_domain_enum = _.unwrap();
        interpreted_program_options.es_momentum_domain = es_momentum_domain_enum_to_variant(es_momentum_domain_enum, raw_program_options.es_n_k);
    } else {
        const std::string message1 = "Problem with interpreting 'es_momentum_domain_string' program option.";
        const std::string message = message1 + "\n" + _.unwrap_err().what();
        throw std::runtime_error(message);
    }
    interpreted_program_options.print_flags.print_unpopulated_basis_flag = raw_program_options.print_unpopulated_basis_flag;
    interpreted_program_options.print_flags.print_unpopulated_basis_size_flag = raw_program_options.print_unpopulated_basis_size_flag;
    interpreted_program_options.print_flags.print_populated_basis_flag = raw_program_options.print_populated_basis_flag;
    interpreted_program_options.print_flags.print_populated_basis_size_flag = raw_program_options.print_populated_basis_size_flag;
    interpreted_program_options.print_flags.print_hamiltonian_stats = raw_program_options.print_hamiltonian_stats;
    interpreted_program_options.print_flags.print_sp_hamiltonian_flag = raw_program_options.print_sp_hamiltonian_flag;
    interpreted_program_options.print_flags.print_hamiltonian_flag = raw_program_options.print_hamiltonian_flag;
    interpreted_program_options.print_flags.print_eigen_values_flag = raw_program_options.print_eigen_values_flag;
    interpreted_program_options.print_flags.print_eigen_vectors_flag = raw_program_options.print_eigen_vectors_flag;
    interpreted_program_options.print_flags.print_pretty_vectors_flag = raw_program_options.print_pretty_vectors_flag;
    interpreted_program_options.print_flags.print_pretty_min_max_n_kstates = {raw_program_options.print_pretty_min_n_kstates, raw_program_options.print_pretty_max_n_kstates};
    interpreted_program_options.print_flags.print_pretty_probability_treshold = raw_program_options.print_pretty_probability_treshold;
    if (raw_program_options.n_threads != 0 && raw_program_options.n_threads <= 256) {
        interpreted_program_options.n_threads = raw_program_options.n_threads;
    } else {
        const std::string message = "Problem with 'n_threads' -- the value must not be equal to 0 and must not be grater than 256.";
        throw std::runtime_error(message);
    }
    return interpreted_program_options;
}
