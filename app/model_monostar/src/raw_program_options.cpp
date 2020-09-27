// METRPOPOLIS:
#include <model_monostar/raw_program_options.hpp>
#include <model_monostar/interpret_model_type_string.hpp>
#include <model_monostar/interpret_run_type_string.hpp>
#include <model_monostar/interpret_es_momentum_domain.hpp>
// EXTENSIONS:
#include <extensions/range_streamer.hpp>
// BOOST:
#include <boost/range/adaptor/map.hpp>
#include <boost/program_options.hpp>
// STD:
#include <iostream>

namespace {

void emit_help(std::ostream& s,
               const boost::program_options::options_description& desc) {
    using namespace extension::boost::stream_pragma;
    s << "Program: model_monostar" << std::endl;
    s << desc << std::endl;
    const auto range_stream_settings = RSS<std::string>().set_null_sustainer().set_string_separer(", ");
    const std::string possible_values_model_type_string = (interpret_model_type_string_map | boost::adaptors::map_keys | range_stream_settings).str();
    const std::string possible_values_run_type_string = (interpret_run_type_string_map | boost::adaptors::map_keys | range_stream_settings).str();
    const std::string possible_values_es_momentum_domain = (interpret_es_momentum_domain_map | boost::adaptors::map_keys | range_stream_settings).str();
    s << "Possible values of model_type string: " << possible_values_model_type_string << "." << std::endl;
    s << "Possible values of run_type string: " << possible_values_run_type_string << "." << std::endl;
    s << "Possible values of es_momentum_domain string: " << possible_values_es_momentum_domain << "." << std::endl;
}

}  // namespace

RawProgramOptions grep_program_options(int argc, char** argv) {
    RawProgramOptions program_options;
    // standard arguments:
    boost::program_options::options_description desc("Options");
    desc.add_options()
            // --help, -h:
            ("help,h", "Print help messages")
            // --n_sites,-n:
            ("n_sites,n",
             boost::program_options::value<unsigned>(&program_options.n_sites)->default_value(8))
            // --n_pt,-p:
            ("n_pt,p",
             boost::program_options::value<unsigned>(&program_options.n_pt)->default_value(2))
            // --model_type_string,-m:
            ("model_type_string,m",
             boost::program_options::value<std::string>(&program_options.model_type_string)->default_value("af"))
            // --hamiltonian_J_classical:
            ("hamiltonian_J_classical",
             boost::program_options::value<double>(&program_options.hamiltonian_J_classical)->default_value(1.0))
            // --hamiltonian_J_quantum:
            ("hamiltonian_J_quantum",
             boost::program_options::value<double>(&program_options.hamiltonian_J_quantum)->default_value(1.0))
            // --hamiltonian_B:
            ("hamiltonian_B",
             boost::program_options::value<double>(&program_options.hamiltonian_B)->default_value(0.0))
            // --hamiltonian_delta:
            ("hamiltonian_delta",
             boost::program_options::value<double>(&program_options.hamiltonian_Pdelta_coef)->default_value(1.0))
            // --hamiltonian_xx:
            ("hamiltonian_xx",
             boost::program_options::value<double>(&program_options.hamiltonian_Pxx_coef)->default_value(1.0))
            // --hamiltonian_xz:
            ("hamiltonian_xz",
             boost::program_options::value<double>(&program_options.hamiltonian_Pxz_coef)->default_value(0.0))
            // --hamiltonian_zz:
            ("hamiltonian_zz",
             boost::program_options::value<double>(&program_options.hamiltonian_Pzz_coef)->default_value(0.0))
            // --theta_opt,o:
            ("theta_opt,o",
             boost::program_options::value<double>(&program_options.theta_opt)->default_value(0.0))
            // --run_type,-r:
            ("run_type,r",
             boost::program_options::value<std::string>(&program_options.run_type_string)->default_value("eg"))
            // --es_momentum_domain_string,-d:
            ("es_momentum_domain_string,d",
             boost::program_options::value<std::string>(&program_options.es_momentum_domain_string)->default_value("all"))
            // --es_n_k,-k:
            ("es_n_k,k",
             boost::program_options::value<unsigned>(&program_options.es_n_k)->default_value(0))
            // --print_unpopulated_basis_flag,-U:
            ("print_unpopulated_basis_flag,U",
             boost::program_options::bool_switch(&program_options.print_unpopulated_basis_flag)->default_value(false))
            // --print_unpopulated_basis_size_flag,-u:
            ("print_unpopulated_basis_size_flag,u",
             boost::program_options::bool_switch(&program_options.print_unpopulated_basis_size_flag)->default_value(false))
            // --print_populated_basis_flag,-B:
            ("print_populated_basis_flag,B",
             boost::program_options::bool_switch(&program_options.print_populated_basis_flag)->default_value(false))
            // --print_populated_basis_size_flag,-b:
            ("print_populated_basis_size_flag,b",
             boost::program_options::bool_switch(&program_options.print_populated_basis_size_flag)->default_value(false))
            // --print_hamiltonian_stats:
            ("print_hamiltonian_stats",
             boost::program_options::bool_switch(&program_options.print_hamiltonian_stats)->default_value(false))
            // --print_sp_hamiltonian_flag,S:
            ("print_sp_hamiltonian_flag,S",
             boost::program_options::bool_switch(&program_options.print_sp_hamiltonian_flag)->default_value(false))
            // --print_hamiltonian_flag,-H:
            ("print_hamiltonian_flag,H",
             boost::program_options::bool_switch(&program_options.print_hamiltonian_flag)->default_value(false))
            // --print_eigen_values_flag,-L:
            ("print_eigen_values_flag,E",
             boost::program_options::bool_switch(&program_options.print_eigen_values_flag)->default_value(false))
            // --print_eigen_vectors_flag,-V:
            ("print_eigen_vectors_flag,V",
             boost::program_options::bool_switch(&program_options.print_eigen_vectors_flag)->default_value(false))
            // --print_pretty_vectors_flag,-i:
            ("print_pretty_vectors_flag,i",
             boost::program_options::bool_switch(&program_options.print_pretty_vectors_flag)->default_value(false))
            // --print_pretty_min_n_kstates:
            ("print_pretty_min_n_kstates",
             boost::program_options::value<unsigned>(&program_options.print_pretty_min_n_kstates)->default_value(5))
            // --print_pretty_max_n_kstates:
            ("print_pretty_max_n_kstates",
             boost::program_options::value<unsigned>(&program_options.print_pretty_max_n_kstates)->default_value(20))
            // --print_pretty_probability_treshold:
            ("print_pretty_probability_treshold",
             boost::program_options::value<double>(&program_options.print_pretty_probability_treshold)->default_value(0.1))
            // --n_threads,j:
            ("n_threads,j",
             boost::program_options::value<unsigned>(&program_options.n_threads)->default_value(4));

    boost::program_options::positional_options_description p;
    // positional arguments:
    p.add("n_sites", 1);
    p.add("n_pt", 1);
    p.add("model", 1);
    boost::program_options::variables_map vm;
    try {
        boost::program_options::store(
                    boost::program_options::command_line_parser(argc, argv).options(desc).positional(p).run(),
                    vm);  // may throw
        if (vm.count("help")) {
            emit_help(std::cout, desc);
            exit(0);
        }
        // sets auto variables (eq. class_specification_file_path),
        // throw is required variable is missing:
        boost::program_options::notify(vm);  // may throw
    } catch (boost::program_options::error& e) {
        std::cerr << "[GLOBAL ERROR] [PROGRAM OPTIONS ERROR]: " << e.what()
                  << std::endl;
        std::cerr << std::endl;
        emit_help(std::cerr, desc);
        exit(1);
    }
    return program_options;
}
