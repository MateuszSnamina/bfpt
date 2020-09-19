#pragma once

#include <model_monostar/enum_es_momentum_domain.hpp>

#include <utility/result.hpp>

#include <map>

#include <stdexcept>
#include <string>

extern const std::map<std::string, EsMomentumDomain> interpret_es_momentum_domain_map;

utility::Result<EsMomentumDomain, std::domain_error> interpret_es_momentum_domain(const std::string&);

EsMomentumDomainVariant es_momentum_domain_enum_to_variant(EsMomentumDomain m, unsigned n_k);
