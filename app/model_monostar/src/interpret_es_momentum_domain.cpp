// SELF:
#include <model_monostar/interpret_es_momentum_domain.hpp>
// EXTENSIONS:
#include <extensions/range_streamer.hpp>
// BOOST:
#include <boost/range/adaptor/map.hpp>

const extern std::map<std::string, EsMomentumDomain> interpret_es_momentum_domain_map{
    {"half_pi_without", EsMomentumDomain::half_pi_without},
    {"half_pi_with", EsMomentumDomain::half_pi_with},
    {"pi_without", EsMomentumDomain::pi_without},
    {"pi_with", EsMomentumDomain::pi_with},
    {"all", EsMomentumDomain::all},
    {"one", EsMomentumDomain::one},
};


utility::Result<EsMomentumDomain, std::domain_error> interpret_es_momentum_domain(const std::string& es_momentum_domain_string) {
    using extension::boost::stream_pragma::RSS;
    using extension::boost::stream_pragma::operator|;
    using extension::boost::stream_pragma::Stringifier;
    using ResultT = utility::Result<EsMomentumDomain, std::domain_error>;
    if (interpret_es_momentum_domain_map.count(es_momentum_domain_string)) {
        return ResultT::Ok(interpret_es_momentum_domain_map.at(es_momentum_domain_string));
    } else {
        const std::string message1 = "Invalid model type string '" + es_momentum_domain_string + "'.";
        const auto range_stream_settings = RSS<std::string>().set_null_sustainer().set_string_separer(", ");
        const std::string possible_values = (interpret_es_momentum_domain_map | boost::adaptors::map_keys | range_stream_settings).str();
        const std::string message2 = "Valid strings are: " + possible_values + ".";
        const std::string message = message1 + " " + message2;
        return ResultT::Err(std::domain_error(message));
    }
}

EsMomentumDomainVariant es_momentum_domain_enum_to_variant(EsMomentumDomain m, unsigned n_k) {
    if (m == EsMomentumDomain::half_pi_without) {
        return EsMomentumDomainHalfPiWithout{};
    } else if (m == EsMomentumDomain::half_pi_with) {
        return EsMomentumDomainHalfPiWith{};
    } else if (m == EsMomentumDomain::pi_without) {
        return EsMomentumDomainPiWithout{};
    } else if (m == EsMomentumDomain::pi_with) {
        return EsMomentumDomainPiWith{};
    } else if (m == EsMomentumDomain::all) {
        return EsMomentumDomainAll{};
    } else if (m == EsMomentumDomain::one) {
        return EsMomentumDomainOne{n_k};
    }  else {
        throw std::logic_error("Invalid model enum value!");
    }
    return EsMomentumDomainAll{};
}
