#pragma once

#include <kstate_impl/kstate_concrete_integral.hpp>

#include <iostream>

using kstate_impl::ctr_from_range;

// #######################################################################
// ## MonostarSiteState -- class                                        ##
// #######################################################################

class MonostarSiteState {
   public:
    constexpr MonostarSiteState(bool is_excited) : _is_excited(is_excited){};
    constexpr bool operator<(const MonostarSiteState& other) const;
    constexpr bool operator==(const MonostarSiteState& other) const;
    constexpr bool operator!=(const MonostarSiteState& other) const;
    friend std::ostream& operator<<(std::ostream&, const MonostarSiteState&);
   //private://TODO RESTORE
    bool _is_excited;
};

// ***********************************************************************

inline constexpr bool MonostarSiteState::operator<(const MonostarSiteState& other) const {
    return static_cast<int>(this->_is_excited) < static_cast<int>(other._is_excited);
}

inline constexpr bool MonostarSiteState::operator==(const MonostarSiteState& other) const {
    return this->_is_excited == other._is_excited;
}

inline constexpr bool MonostarSiteState::operator!=(const MonostarSiteState& other) const {
    return this->_is_excited != other._is_excited;
}

inline std::ostream&
operator<<(std::ostream& stream, const MonostarSiteState& state) {
    stream << (state._is_excited ? '*' : '_');
    return stream;
}

// #######################################################################
// ## MonostarSiteState -- global values                                ##
// #######################################################################

//TODO constexpr
constexpr MonostarSiteState gs{false};  // ground-state (there is no star)
constexpr MonostarSiteState es{true};  // excited-state (there is a star)

// #######################################################################
// ## MonostarSiteState -- implement trait                              ##
// #######################################################################

namespace kstate_trait {

template<>
struct TraitSiteState<MonostarSiteState> {
    static constexpr bool is_site_state_trait = true;
    using SiteStateT = MonostarSiteState;

    constexpr static unsigned site_basis_dim() {
        return 2u;
    }

    constexpr static unsigned get_index(const SiteStateT& state) {
        if (state == gs) {
            return 0u;
        } else if (state == es) {
            return 1u;
        } else {
            throw std::domain_error("Not a valid state."); //TODO restore
        }
    }

    constexpr static SiteStateT from_index(unsigned idx) {
        if (idx == 0u) {
            return gs;
        } else if (idx == 1u) {
            return es;
        } else {
            throw std::domain_error("Index out of range.");
        }
    }
};

} // end of namespace kstate


using MonostarSiteStateTrait = kstate_trait::TraitSiteState<MonostarSiteState>;
