#pragma once

#include <kstate_trait/trait_site_state.hpp>

#include <boost/operators.hpp>

#include <iostream>
#include <vector>

// #######################################################################
// ## MonostarSiteState -- class                                        ##
// #######################################################################

class MySiteState : boost::totally_ordered<MySiteState> {
   public:
    constexpr MySiteState(int i) : _i(i){};
    constexpr bool operator<(const MySiteState& other) const;
    constexpr bool operator==(const MySiteState& other) const;
    friend std::ostream& operator<<(std::ostream&, const MySiteState&);
    int _i;
};

// ***********************************************************************

inline constexpr bool MySiteState::operator<(const MySiteState& other) const {
    return this->_i < other._i;
}

inline constexpr bool MySiteState::operator==(const MySiteState& other) const {
    return this->_i == other._i;
}

inline std::ostream&
operator<<(std::ostream& stream, const MySiteState& state) {
    stream << ('A' + state._i);
    return stream;
}

// #######################################################################
// ## MonostarSiteState -- implement trait                              ##
// #######################################################################

namespace kstate_trait {

template<>
struct TraitSiteState<MySiteState> {
    static constexpr bool is_site_state_trait = true;
    using SiteStateT = MySiteState;

    constexpr static unsigned site_basis_dim() {
        return 5u;
    }

    constexpr static unsigned get_index(const SiteStateT& state) {
        if (state._i >= 10 && state._i < 15) {
            return state._i - 10;
        } else {
            throw std::domain_error("Not a valid state."); //TODO restore
        }
    }

    constexpr static SiteStateT from_index(unsigned idx) {
        if (idx < 5) {
            return SiteStateT(idx + 10);
        } else {
            throw std::domain_error("Index out of range.");
        }
    }
};

} // end of namespace kstate

using MySiteStateTrait = kstate_trait::TraitSiteState<MySiteState>;

