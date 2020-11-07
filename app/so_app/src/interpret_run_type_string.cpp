// SELF:
#include <so_app/interpret_run_type_string.hpp>
// EXTENSIONS:
#include <extensions/range_streamer.hpp>
// BOOST:
#include <boost/range/adaptor/map.hpp>

namespace so_app {

const extern std::map<std::string, RunType> interpret_run_type_string_map{
    {"g", RunType::G},
    {"e", RunType::E},
    {"eg", RunType::EG},
};

utility::Result<RunType, std::domain_error> interpret_run_type_string(const std::string& run_type_string) {
    using extension::boost::stream_pragma::RSS;
    using extension::boost::stream_pragma::operator|;
    using extension::boost::stream_pragma::Stringifier;
    using ResultT = utility::Result<RunType, std::domain_error>;
    if (interpret_run_type_string_map.count(run_type_string)) {
        return ResultT::Ok(interpret_run_type_string_map.at(run_type_string));
    } else {
        const std::string message1 = "Invalid run type string '" + run_type_string + "'.";
        const auto range_stream_settings = RSS<std::string>().set_null_sustainer().set_string_separer(", ");
        const std::string possible_values = (interpret_run_type_string_map | boost::adaptors::map_keys | range_stream_settings).str();
        const std::string message2 = "Valid strings are: " + possible_values + ".";
        const std::string message = message1 + " " + message2;
        return ResultT::Err(std::domain_error(message));
    }
}

}  // end of namespace so_app
