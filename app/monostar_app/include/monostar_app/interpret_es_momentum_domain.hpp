#pragma once

#include <monostar_app/enum_es_momentum_domain.hpp>

#include <utility/result.hpp>

#include <map>

#include <stdexcept>
#include <string>

namespace monostar_app {

extern const std::map<std::string, EsMomentumDomain> interpret_es_momentum_domain_map;

utility::Result<EsMomentumDomain, std::domain_error> interpret_es_momentum_domain(const std::string&);

EsMomentumDomainVariant es_momentum_domain_enum_to_variant(EsMomentumDomain m, unsigned n_k);

} // end of namespace monostar_app
