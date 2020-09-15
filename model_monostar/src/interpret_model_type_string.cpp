// SELF:
#include <model_monostar/interpret_model_type_string.hpp>
// EXTENSIONS:
#include <extensions/range_streamer.hpp>
// BOOST:
#include <boost/range/adaptor/map.hpp>

const extern std::map<std::string, ModelType> interpret_model_type_string_map{
    {"af", ModelType::AF},
    {"fm", ModelType::FM},
};

utility::Result<ModelType, std::domain_error> interpret_model_type_string(const std::string& model_type_string) {
    using extension::boost::stream_pragma::RSS;
    using extension::boost::stream_pragma::operator|;
    using extension::boost::stream_pragma::Stringifier;
    //using namespace extension::boost::stream_pragma;//  TEMP TODO remove
    using ResultT = utility::Result<ModelType, std::domain_error>;
    if (interpret_model_type_string_map.count(model_type_string)) {
        return ResultT::Ok(interpret_model_type_string_map.at(model_type_string));
    } else {
        const std::string message1 = "Invalid model type string '" + model_type_string + "'.";
        const auto range_stream_settings = RSS<std::string>().set_null_sustainer().set_string_separer(", ");
        const std::string possible_values = (interpret_model_type_string_map | boost::adaptors::map_keys | range_stream_settings).str();
        const std::string message2 = "Valid strings are: " + possible_values + ".";
        const std::string message = message1 + " " + message2;
        return ResultT::Err(std::domain_error(message));
    }
}
