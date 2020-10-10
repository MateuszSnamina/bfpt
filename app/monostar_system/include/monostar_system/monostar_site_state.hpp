#pragma once

#include <kstate/trait_site_state.hpp>

#include <boost/operators.hpp>

#include <iostream>
#include <vector>

// #######################################################################
// ## MonostarSiteState -- class                                        ##
// #######################################################################

namespace monostar_system {

class MonostarSiteState : boost::totally_ordered<MonostarSiteState> {
   public:
    MonostarSiteState(bool is_excited) : _is_excited(is_excited){};
    bool operator<(const MonostarSiteState& other) const;
    bool operator==(const MonostarSiteState& other) const;
    friend std::ostream& operator<<(std::ostream&, const MonostarSiteState&);
   //private://RESTORE
    bool _is_excited;
};

// ***********************************************************************

inline bool MonostarSiteState::operator<(const MonostarSiteState& other) const {
    return static_cast<int>(this->_is_excited) < static_cast<int>(other._is_excited);
}

inline bool MonostarSiteState::operator==(const MonostarSiteState& other) const {
    return this->_is_excited == other._is_excited;
}

inline std::ostream&
operator<<(std::ostream& stream, const MonostarSiteState& state) {
    stream << (state._is_excited ? '*' : '_');
    return stream;
}

}  // namespace monostar_system

// #######################################################################
// ## MonostarSiteState -- global values                                ##
// #######################################################################

namespace monostar_system {

extern const MonostarSiteState gs;  // ground-state (there is no star)
extern const MonostarSiteState es;  // excited-state (there is a star)
extern const std::vector<MonostarSiteState> ordered_site_states;//TODO remove

}  // namespace monostar_system

// #######################################################################
// ## MonostarSiteState -- implement trait                              ##
// #######################################################################

namespace kstate {

template<>
struct TraitSiteState<monostar_system::MonostarSiteState> {
    static constexpr bool is_site_state_trait = true;
    using SiteStateT = monostar_system::MonostarSiteState;

    static unsigned site_basis_dim() {
        return 2u;
    }

    static unsigned get_index(const SiteStateT& state) {
        if (state == monostar_system::gs) {
            return 0u;
        } else if (state == monostar_system::es) {
            return 1u;
        } else {
            //throw std::domain_error("Not a valid state."); //TODO restore
            throw std::domain_error("Not a valid state.NNNN"); //TODO remove
        }
    }

    static SiteStateT from_index(unsigned idx) {
        if (idx == 0u) {
            return monostar_system::gs;
        } else if (idx == 1u) {
            return monostar_system::es;
        } else {
            throw std::domain_error("Index out of range.");
        }
    }
};

} // end of namespace kstate

namespace monostar_system {

using MonostarSiteStateTrait = kstate::TraitSiteState<monostar_system::MonostarSiteState>;

} // end of namespace monostar_system
