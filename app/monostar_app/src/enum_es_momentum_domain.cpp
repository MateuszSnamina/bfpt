#include<monostar_app/enum_es_momentum_domain.hpp>

#include <stdexcept>

// -----------------------------------------------------
// --  Helpers                                        --
// -----------------------------------------------------

namespace {

using namespace monostar_app;

class EsMomentumDomainVariantStreamerVisitor {
public:
    EsMomentumDomainVariantStreamerVisitor(std::ostream& stream)
        : _stream(stream) {
    }

    void operator() (const EsMomentumDomainHalfPiWithout&) {
        _stream << "EsMomentumDomainHalfPiWithout";
    }

    void operator() (const EsMomentumDomainHalfPiWith&) {
        _stream << "EsMomentumDomainHalfPiWith";
    }

    void operator() (const EsMomentumDomainPiWithout&) {
        _stream << "EsMomentumDomainPiWithout";
    }

    void operator() (const EsMomentumDomainPiWith&) {
        _stream << "EsMomentumDomainPiWith";
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


class EsMomentumDomainVariantToMomentumRangeSapnVisitor {
public:
    EsMomentumDomainVariantToMomentumRangeSapnVisitor(unsigned n_sites)
        : _n_sites(n_sites) {
    }

    std::pair<unsigned, unsigned> operator() (const EsMomentumDomainHalfPiWithout&) const {
        return {0, (_n_sites + 3) / 4};
    }

    std::pair<unsigned, unsigned> operator() (const EsMomentumDomainHalfPiWith&) const  {
        return {0, (_n_sites + 4) / 4};
    }

    std::pair<unsigned, unsigned> operator() (const EsMomentumDomainPiWithout&) const {
        return {0, (_n_sites + 1) / 2};
    }

    std::pair<unsigned, unsigned> operator() (const EsMomentumDomainPiWith&) const {
        return {0, (_n_sites + 2) / 2};
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

}

// -----------------------------------------------------
// --  Enum type                                      --
// -----------------------------------------------------

namespace monostar_app {

std::ostream& operator<<(std::ostream& s, const EsMomentumDomain& m) {
    if (m == EsMomentumDomain::half_pi_without) {
        s << "EsMomentumDomain::half_pi_without";
    } else if (m == EsMomentumDomain::half_pi_with) {
        s << "EsMomentumDomain::half_pi_with";
    } else if (m == EsMomentumDomain::pi_without) {
        s << "EsMomentumDomain::pi_without";
    } else if (m == EsMomentumDomain::pi_with) {
        s << "EsMomentumDomain::pi_with";
    } else if (m == EsMomentumDomain::all) {
        s << "EsMomentumDomain::all";
    } else if (m == EsMomentumDomain::one) {
        s << "EsMomentumDomain::one";
    } else {
        throw std::logic_error("Invalid model enum value!");
    }
    return s;
}

} // end of namespace monostar_app

// -----------------------------------------------------
// --  Variant type                                   --
// -----------------------------------------------------

namespace monostar_app {

std::ostream&
operator<<(std::ostream& s, const EsMomentumDomainVariant& v) {
    EsMomentumDomainVariantStreamerVisitor visitor{s};
    std::visit(visitor, v);
    return s;
}

std::pair<unsigned, unsigned>
es_momentum_domain_variant_to_momentum_range_sapn(EsMomentumDomainVariant md, unsigned n_sites) {
    const EsMomentumDomainVariantToMomentumRangeSapnVisitor visitor{n_sites};
    return std::visit(visitor, md);
}

} // end of namespace monostar_app
