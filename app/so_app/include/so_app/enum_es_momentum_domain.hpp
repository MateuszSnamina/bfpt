#pragma once

// STD:
#include <ostream>
#include <variant>

namespace so_app {

// -----------------------------------------------------
// --  Enum type                                      --
// -----------------------------------------------------

enum class EsMomentumDomain {
    half_pi_without,
    half_pi_with,
    pi_without,
    pi_with,
    all,
    one
};

std::ostream& operator<<(std::ostream&, const EsMomentumDomain&);

// -----------------------------------------------------
// --  Variant type                                   --
// -----------------------------------------------------

struct EsMomentumDomainHalfPiWithout {
};

struct EsMomentumDomainHalfPiWith {
};

struct EsMomentumDomainPiWithout {
};

struct EsMomentumDomainPiWith {
};

struct EsMomentumDomainAll {
};

struct EsMomentumDomainOne {
    unsigned n_k;
};

using EsMomentumDomainVariant =
    std::variant<
        EsMomentumDomainHalfPiWithout,
        EsMomentumDomainHalfPiWith,
        EsMomentumDomainPiWithout,
        EsMomentumDomainPiWith,
        EsMomentumDomainAll,
        EsMomentumDomainOne>;

std::ostream&
operator<<(std::ostream&, const EsMomentumDomainVariant&);

std::pair<unsigned, unsigned>
es_momentum_domain_variant_to_momentum_range_sapn(EsMomentumDomainVariant md, unsigned n_sites);

}  // end of namespace so_app
