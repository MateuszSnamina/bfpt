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
             boost::program_options::value<double>(&program_options.J_quantum)->default_value(1.0));
    boost::program_options::positional_options_description p;
    // positional arguments:
    p.add("temperature_steps", 1);
    p.add("n_thermal", 1);
    p.add("n_average", 1);
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
