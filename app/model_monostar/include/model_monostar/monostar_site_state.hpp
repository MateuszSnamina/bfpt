#ifndef MODEL_MONOSTAR_MONOSTAR_SITE_STATE_HPP
#define MODEL_MONOSTAR_MONOSTAR_SITE_STATE_HPP

#include <kstate/trait_site_state.hpp>

#include <boost/operators.hpp>

#include <iostream>
#include <vector>

// #######################################################################
// ## MonostarSiteState -- class                                        ##
// #######################################################################

namespace model_monostar {

class MonostarSiteState : boost::totally_ordered<MonostarSiteState> {
   public:
    MonostarSiteState(bool is_excited) : _is_excited(is_excited){};
    bool operator<(const MonostarSiteState& other) const;
    bool operator==(const MonostarSiteState& other) const;
    friend std::ostream& operator<<(std::ostream&, const MonostarSiteState&);
   private:
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

}  // namespace model_monostar

// #######################################################################
// ## MonostarSiteState -- global values                                ##
// #######################################################################

namespace model_monostar {

extern const MonostarSiteState gs;  // ground-state (there is no star)
extern const MonostarSiteState es;  // excited-state (there is a star)
extern const std::vector<MonostarSiteState> ordered_site_states;

}  // namespace model_monostar

// #######################################################################
// ## MonostarSiteState -- implement trait                              ##
// #######################################################################

namespace kstate {

template<>
struct TraitSiteState<model_monostar::MonostarSiteState> {
    static constexpr bool is_site_state_trait = true;
    using SiteStateT = model_monostar::MonostarSiteState;

    static unsigned site_basis_dim() {
        return 2u;
    }

    static unsigned get_index(const SiteStateT& state) {
        if (state == model_monostar::gs) {
            return 0u;
        } else if (state == model_monostar::es) {
            return 1u;
        } else {
            throw std::domain_error("Not a valid state.");
        }
    }

    static SiteStateT from_index(unsigned idx) {
        if (idx == 0u) {
            return model_monostar::gs;
        } else if (idx == 1u) {
            return model_monostar::es;
        } else {
            throw std::domain_error("Index out of range.");
        }
    }
};

} // end of namespace kstate

namespace model_monostar {

using MonostarSiteStateTrait = kstate::TraitSiteState<model_monostar::MonostarSiteState>;

} // end of namespace model_monostar



#endif
