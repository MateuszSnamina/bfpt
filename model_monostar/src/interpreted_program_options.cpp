#include <model_monostar/interpreted_program_options.hpp>

#include <model_monostar/raw_program_options.hpp>
#include <model_monostar/interpret_model_type_string.hpp>

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
    interpreted_program_options.J_classical = raw_program_options.J_classical;
    interpreted_program_options.J_quantum = raw_program_options.J_quantum;
    return interpreted_program_options;
}
