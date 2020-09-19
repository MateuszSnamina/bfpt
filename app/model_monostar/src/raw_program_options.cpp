// METRPOPOLIS:
#include <model_monostar/raw_program_options.hpp>
#include <model_monostar/interpret_model_type_string.hpp>
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
    const std::string possible_values_model_type_string =
            (interpret_model_type_string_map | boost::adaptors::map_keys | range_stream_settings).str();
    std::cout << "Possible values of model_type string: " << possible_values_model_type_string << "." << std::endl;
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
            // --J_classical,-c:
            ("J_classical,c",
             boost::program_options::value<double>(&program_options.J_classical)->default_value(1.0))
            // --J_quantum,-q:
            ("J_quantum,q",
             boost::program_options::value<double>(&program_options.J_quantum)->default_value(1.0))
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
            // --print_unpopulated_basis_size_flag:
            ("print_unpopulated_basis_size_flag",
             boost::program_options::bool_switch(&program_options.print_unpopulated_basis_size_flag)->default_value(true))
            // --print_populated_basis_flag,-b:
            ("print_populated_basis_flag,b",
             boost::program_options::bool_switch(&program_options.print_populated_basis_flag)->default_value(false))
            // --print_populated_basis_size_flag:
            ("print_populated_basis_size_flag",
             boost::program_options::bool_switch(&program_options.print_populated_basis_size_flag)->default_value(true))
            // --print_sp_hamiltonian_flag,S:
            ("print_sp_hamiltonian_flag,S",
             boost::program_options::bool_switch(&program_options.print_sp_hamiltonian_flag)->default_value(false))
            // --print_hamiltonian_flag,-H:
            ("print_hamiltonian_flag,H",
             boost::program_options::bool_switch(&program_options.print_hamiltonian_flag)->default_value(false))
            // --print_eigen_values_flag,-L:
            ("print_eigen_values_flag,E",
             boost::program_options::bool_switch(&program_options.print_eigen_values_flag)->default_value(true))
            // --print_eigen_vectors_flag,-V:
            ("print_eigen_vectors_flag,V",
             boost::program_options::bool_switch(&program_options.print_eigen_vectors_flag)->default_value(false));
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
