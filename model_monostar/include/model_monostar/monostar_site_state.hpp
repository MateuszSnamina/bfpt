#ifndef MODEL_MONOSTAR_MONOSTAR_SITE_STATE_HPP
#define MODEL_MONOSTAR_MONOSTAR_SITE_STATE_HPP

#include <model_monostar/monostar_site_state.hpp>

#include <boost/operators.hpp>
#include <iostream>

// #######################################################################
// ## MonostarSiteState -- class                                        ##
// #######################################################################

namespace model_monostar {

class MonostarSiteState : boost::totally_ordered<MonostarSiteState> {
   public:
    MonostarSiteState(bool is_excited) : _is_excited(is_excited){};
    bool operator<(const MonostarSiteState& other) const;
    bool operator==(const MonostarSiteState& other) const;

   public:  // TODO: make it private
    bool _is_excited;
};

// ***********************************************************************

inline bool MonostarSiteState::operator<(const MonostarSiteState& other) const {
    return static_cast<int>(this->_is_excited) < static_cast<int>(other._is_excited);
}

inline bool MonostarSiteState::operator==(const MonostarSiteState& other) const {
    return static_cast<int>(this->_is_excited) == static_cast<int>(other._is_excited);
}

}  // namespace model_monostar

// #######################################################################
// ## MonostarSiteState -- global values                                ##
// #######################################################################

namespace model_monostar {

extern const MonostarSiteState gs;  // ground-state (there is no star)
extern const MonostarSiteState es;  // excited-state (there is a star)

}  // namespace model_monostar

// #######################################################################
// ## MonostarSiteState -- printing                                     ##
// #######################################################################

namespace model_monostar {

inline std::ostream&
operator<<(std::ostream& stream, const MonostarSiteState& state) {
    stream << (state._is_excited ? '*' : '_');
    return stream;
}

}  // namespace model_monostar

#endif