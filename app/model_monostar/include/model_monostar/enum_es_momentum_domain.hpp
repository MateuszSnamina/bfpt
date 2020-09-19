#pragma once

// STD:
#include <string>
#include <ostream>
#include <stdexcept>
#include <variant>

// -----------------------------------------------------
// --  Enum type                                      --
// -----------------------------------------------------

enum class EsMomentumDomain {
    quarter_pi_without,
    quarter_pi_with,
    half_pi_without,
    half_pi_with,
    all,
    one
};

inline std::ostream& operator<<(std::ostream& s, const EsMomentumDomain& m) {
    if (m == EsMomentumDomain::quarter_pi_without) {
        s << "EsMomentumDomain::quarter_pi_without";
    } else if (m == EsMomentumDomain::quarter_pi_with) {
        s << "EsMomentumDomain::quarter_pi_with";
    } else if (m == EsMomentumDomain::half_pi_without) {
        s << "EsMomentumDomain::half_pi_without";
    } else if (m == EsMomentumDomain::half_pi_with) {
        s << "EsMomentumDomain::half_pi_with";
    } else if (m == EsMomentumDomain::all) {
        s << "EsMomentumDomain::all";
    } else if (m == EsMomentumDomain::one) {
        s << "EsMomentumDomain::one";
    }  else {
        throw std::logic_error("Invalid model enum value!");
    }
    return s;
}

// -----------------------------------------------------
// --  Variant type                                   --
// -----------------------------------------------------


struct EsMomentumDomainQuarterPiWithout {
};

struct EsMomentumDomainQuarterPiWith {
};

struct EsMomentumDomainHalfPiWithout {
};

struct EsMomentumDomainHalfPiWith {
};

struct EsMomentumDomainAll {
};

struct EsMomentumDomainOne {
    unsigned n_k;
};

using EsMomentumDomainVariant =
std::variant<
EsMomentumDomainQuarterPiWithout,
EsMomentumDomainQuarterPiWith,
EsMomentumDomainHalfPiWithout,
EsMomentumDomainHalfPiWith,
EsMomentumDomainAll,
EsMomentumDomainOne
>;

class EsMomentumDomainVariantStreamerVisitor {
public:
    EsMomentumDomainVariantStreamerVisitor(std::ostream& stream)
        : _stream(stream) {
    }

    void operator() (const EsMomentumDomainQuarterPiWithout&) {
        _stream << "EsMomentumDomainQuarterPiWithout";
    }

    void operator() (const EsMomentumDomainQuarterPiWith&) {
        _stream << "EsMomentumDomainQuarterPiWith";
    }

    void operator() (const EsMomentumDomainHalfPiWithout&) {
        _stream << "EsMomentumDomainHalfPiWithout";
    }

    void operator() (const EsMomentumDomainHalfPiWith&) {
        _stream << "EsMomentumDomainHalfPiWith";
    }

    void operator() (const EsMomentumDomainAll&) {
        _stream << "EsMomentumDomainAll";
    }

    void operator() (const EsMomentumDomainOne& o) {
        _stream << "EsMomentumDomainOne(n_k:" << o.n_k << ")";
    }

private:
    std::ostream& _stream;
};

inline
std::ostream&
operator<<(std::ostream& s, const EsMomentumDomainVariant& v) {
    EsMomentumDomainVariantStreamerVisitor visitor{s};
    std::visit(visitor, v);
    return s;
}

class EsMomentumDomainVariantToMomentumRangeSapnVisitor {
public:
    EsMomentumDomainVariantToMomentumRangeSapnVisitor(unsigned n_sites)
        : _n_sites(n_sites) {
    }

    std::pair<unsigned, unsigned> operator() (const EsMomentumDomainQuarterPiWithout&) const {
        return {0, _n_sites / 4};
    }

    std::pair<unsigned, unsigned> operator() (const EsMomentumDomainQuarterPiWith&) const  {
        return {0, (_n_sites + 1) / 4};
    }

    std::pair<unsigned, unsigned> operator() (const EsMomentumDomainHalfPiWithout&) const {
        return {0, _n_sites / 2};
    }

    std::pair<unsigned, unsigned> operator() (const EsMomentumDomainHalfPiWith&) const {
        return {0, (_n_sites + 1) / 2};
    }

    std::pair<unsigned, unsigned> operator() (const EsMomentumDomainAll&) const {
        return {0, _n_sites};
    }

    std::pair<unsigned, unsigned> operator() (const EsMomentumDomainOne& o) const {
        return {o.n_k, o.n_k + 1};
    }

private:
    unsigned _n_sites;
};

inline
std::pair<unsigned, unsigned>
es_momentum_domain_variant_to_momentum_range_sapn(EsMomentumDomainVariant md, unsigned n_sites) {
    const EsMomentumDomainVariantToMomentumRangeSapnVisitor visitor{n_sites};
    return std::visit(visitor, md);
}
